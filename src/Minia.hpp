/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _GATB_TOOLS_MINIA_HPP_
#define _GATB_TOOLS_MINIA_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

/********************************************************************************/

/**
 */
class Minia : public gatb::core::tools::misc::impl::Tool
{
public:

    /** Constructor. */
    Minia ();

private:

    /** \copydoc Tool::execute. */
    void  execute ();

    /** */
    void assemble (const Graph& graph);
    
    /** */
    void assembleFrom (Node startingNode, Traversal *traversal, const Graph& graph, IBank *outputBank);

    /** */
    void buildSequence (
        const Graph& graph,
        const Node& startingNode,
        size_t length,
        size_t nbContigs,
        const Path& consensusRight,
        const Path& consensusLeft,
        Sequence& seq
    );

    bool keepIsolatedTigs;
    u_int64_t nbContigs         ;
    u_int64_t nbSmallContigs    ;
    u_int64_t totalNt           ;
    u_int64_t maxContigLen      ;
    u_int64_t maxContigLenLeft  ;
    u_int64_t maxContigLenRight ;

};

/********************************************************************************/

#endif /* _GATB_TOOLS_MINIA_HPP_ */

