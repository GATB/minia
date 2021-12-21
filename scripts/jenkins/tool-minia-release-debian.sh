#!/bin/bash
#--------------------------------------------------------------#
#         Continuous integration script for Jenkins            #
#--------------------------------------------------------------#
#
# Default mode :
# This script will exit with error (exit code 1) if any of its steps fails.
# To change this behaviour, choose DO_NOT_STOP_AT_ERROR in Jenkins (see below).
#--------------------------------------------------------------#
set +xv

echo "
-----------------------------------------
 Miscellaneous information 
-----------------------------------------
date      : `date`
hostname  : `hostname`
pwd       : `pwd`

-----------------------------------------
 Jenkins build parameters (user defined)
-----------------------------------------
BRANCH_TO_BUILD      : ${BRANCH_TO_BUILD}
RELEASE_TO_BUILD     : ${RELEASE_TO_BUILD}
INRIA_FORGE_LOGIN    : ${INRIA_FORGE_LOGIN}
TEST_VARIABLE        : ${TEST_VARIABLE}
DO_NOT_STOP_AT_ERROR : ${DO_NOT_STOP_AT_ERROR}

-----------------------------------------
 Jenkins build parameters (built in)
-----------------------------------------
BUILD_NUMBER         : ${BUILD_NUMBER}
JENKINS_HOME         : ${JENKINS_HOME}
WORKSPACE            : ${WORKSPACE}
"
set -xv

# quick look at resources
#-----------------------------------------------
free -h
#-----------------------------------------------
lstopo
#-----------------------------------------------
df -kh
#-----------------------------------------------


################################################################
#                       PREPARE RELEASE                        #
################################################################

# paths to access tool source code and build
JENKINS_TASK=tool-${TOOL_NAME}-build-debian7-64bits-gcc-4.7
BUILD_DIR=/scratchdir/$JENKINS_TASK/gatb-${TOOL_NAME}-release
TOOL_GIT_HOME="/scratchdir/builds/workspace/gatb-${TOOL_NAME}"

# path to 'github_release_manager.sh' script
GRM_PATH="${BUILD_DIR}/github-release-api"
GRM_CMD="${GRM_PATH}/github_release_manager.sh"
# github credentials and repository
GITHUB_REPO=${TOOL_NAME}
GITHUB_OWNER=GATB
GRM_CREDENTIALS="-l $GITHUB_ADMIN -t $GITHUB_TOKEN -o ${GITHUB_OWNER} -r ${GITHUB_REPO}"

# Prepare build dir
rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR

#-----------------------------------------------
# check tag version; 'master' is not allowed
if [ ! "${BRANCH_TO_BUILD}" == "master" ] ; then
	cd ${TOOL_GIT_HOME}
    DOES_TAG_EXIST=`git tag -l | grep "^${BRANCH_TO_BUILD}$"`
    if [ -z ${DOES_TAG_EXIST} ] ; then
    	echo "/!\ Error: tag '${BRANCH_TO_BUILD}' does not exist on 'gatb-${TOOL_NAME}' repository"
        exit 1
    fi
else
    echo "/!\ Error: cannot make an official release on 'master' branch"
    exit 1
fi

#-----------------------------------------------
if [ "$INRIA_FORGE_LOGIN" == none ]; then
	echo "/!\ Error: No login name to connect to Inria Forge"
    exit 1
fi

cd $BUILD_DIR
git clone https://github.com/GATB/github-release-api.git

################################################################
#                       RETRIEVE ARCHIVES FROM INRIA FORGE     #
################################################################

CI_URL=https://ci.inria.fr/gatb-core/view/Minia-gitlab/job
JENKINS_TASK_DEB=tool-minia-build-debian7-64bits-gcc-4.7-gitlab
JENKINS_TASK_MAC=tool-minia-build-macos-10.9.5-gcc-4.2.1-gitlab

#retrieve last build from ci-inria (see tool-lean-build-XXX tasks)
wget $CI_URL/$JENKINS_TASK_DEB/lastSuccessfulBuild/artifact/$JENKINS_TASK_DEB/${TOOL_NAME}-${BRANCH_TO_BUILD}-bin-Linux.tar.gz
[ $? != 0 ] && exit 1

#wget $CI_URL/$JENKINS_TASK_MAC/lastSuccessfulBuild/artifact/${TOOL_NAME}-${BRANCH_TO_BUILD}-bin-Darwin.tar.gz
#[ $? != 0 ] && exit 1 # disabled because mavericks machine is down in 2021

wget $CI_URL/$JENKINS_TASK_DEB/lastSuccessfulBuild/artifact/$JENKINS_TASK_DEB/${TOOL_NAME}-${BRANCH_TO_BUILD}-Source.tar.gz
[ $? != 0 ] && exit 1

################################################################
#                       INTERACT WITH GITHUB                   #
################################################################

# create Github release
${GRM_CMD} ${GRM_CREDENTIALS} -d ${BRANCH_TO_BUILD} -c create
if [ $? != 0 ] ; then
  echo "/!\ Error: unable to create release, check above error"
  exit 1
fi

#upload files
function uploadFile(){
  local FILE_TO_LOAD=$1
  echo "Uploading: ${FILE_TO_LOAD}"
  ${GRM_CMD} ${GRM_CREDENTIALS} -d ${BRANCH_TO_BUILD} -c upload ${FILE_TO_LOAD}
  if [ $? != 0 ] ; then
    echo "/!\ Error: unable to upload file, check above error"
    exit 1
  fi
}

uploadFile ${TOOL_NAME}-${BRANCH_TO_BUILD}-bin-Linux.tar.gz
#uploadFile ${TOOL_NAME}-${BRANCH_TO_BUILD}-bin-Darwin.tar.gz
uploadFile ${TOOL_NAME}-${BRANCH_TO_BUILD}-Source.tar.gz


