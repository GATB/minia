/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <Utils.hpp>

using namespace std;

/********************************************************************************/

#define DEBUG(a)   //a

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
float needleman_wunch (const Path& a, const Path& b)
{
    float gap_score = -5;
    float mismatch_score = -5;
    float match_score = 10;
    #define nw_score(x,y) ( (x == y) ? match_score : mismatch_score )

    int n_a = a.size(), n_b = b.size();
    float ** score =  (float **) malloc (sizeof(float*) * (n_a+1));
    for (int ii=0; ii<(n_a+1); ii++)
    {
        score [ii] = (float *) malloc (sizeof(float) * (n_b+1));
    }

   // float score[n_a+1][n_b+1];  //stack is too small
    // float pointer[n_a+1][n_b+1];

    for (int i = 0; i <= n_a; i++)
        score[i][0] = gap_score * i;
    for (int j = 0; j <= n_b; j++)
        score[0][j] = gap_score * j;

    // compute dp
    for (int i = 1; i <= n_a; i++)
    {
        for (int j = 1; j <= n_b; j++)
        {
            float match = score[i - 1][j - 1] + nw_score(a[i-1],b[j-1]);
            float del =  score[i - 1][j] + gap_score;
            float insert = score[i][j - 1] + gap_score;
            score[i][j] = max( max(match, del), insert);
        }
    }

    // traceback
    int i=n_a, j=n_b;
    float identity = 0;
    while (i > 0 && j > 0)
    {
        float score_current = score[i][j], score_diagonal = score[i-1][j-1], score_up = score[i][j-1], score_left = score[i-1][j];
        if (score_current == score_diagonal + nw_score(a[i-1], b[j-1]))
        {
            if (a[i-1]== b[j-1])
                identity++;
            i -= 1;
            j -= 1;
        }
        else
        {
            if (score_current == score_left + gap_score)
                i -= 1;
            else if (score_current == score_up + gap_score)
                    j -= 1;
        }
    }
    identity /= max( n_a, n_b); // modif GR 27/09/2013    max of two sizes, otherwise free gaps



    for (int ii=0; ii<(n_a+1); ii++)
    {
        free (score [ii]);
    }
    free(score);

    return identity;
}
