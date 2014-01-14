/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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

#ifndef _GATB_TOOLS_TERMINATOR_HPP_
#define _GATB_TOOLS_TERMINATOR_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

/********************************************************************************/

template <typename Item> class ListSet
{
public:
    std::vector<Item> liste;
public:
    void insert (const Item& elem)  {     liste.push_back(elem); }
    void finalize() {    sort(liste.begin(), liste.end()); }
    bool contains(const Item& elem)  const {    return binary_search(liste.begin(), liste.end(), elem);  }
    uint64_t capacity() {return (u_int64_t)liste.capacity();}
    static const int bits_per_element = sizeof(Item)*8;
    ListSet(u_int64_t taille_approx)  {   liste.reserve(taille_approx);  }
    ListSet()  {}
};

/********************************************************************************/

template <typename Key, typename Value>  class AssocSet : public ListSet<Key>
{
public:
    std::vector<Value> liste_value;

public:
    int get (const Key& elem, Value& val) const
    {
        typename std::vector<Key>::const_iterator it;
        it = lower_bound(this->liste.begin(), this->liste.end(),elem);
        if (it == this->liste.end() || elem != *it) return 0;
        size_t rank = it - this->liste.begin();
        val = liste_value[rank];
        return 1;
    }

    int set (const Key& elem, const Value& val)
    {
        typename  std::vector<Key>::iterator it;
        it = lower_bound(this->liste.begin(), this->liste.end(),elem);
        if (it == this->liste.end() ||elem != *it) return 0;
        size_t rank = it - this->liste.begin();
        liste_value[rank]=val;
        return 1;
    }

    void finalize (bool doSort=true)
    {
        if (doSort) {  sort(this->liste.begin(), this->liste.end());  }
        liste_value.assign(this->liste.size(),0);
    }

    AssocSet() {}

    void clear()  { liste_value.assign(liste_value.size(),0); }


    void start_iterator()  {   iterator = this->liste.begin()-1;  }
    bool next_iterator()
    {
        iterator++;
        if (iterator==this->liste.end()) return false;
        return true;
    }

    typename std::vector<Key>::iterator iterator;
};

/********************************************************************************/

class Terminator
{
public:
    Terminator (const Graph& graph) : _graph(graph)  {}

    virtual ~Terminator ()  {}

    const Graph& getGraph() const { return _graph; }

    virtual bool isEnabled () const { return true; }

    virtual void mark      (const Edge& edge) = 0;
    virtual bool is_marked (const Edge& edge)  const = 0;

    virtual void mark      (const Node& node) = 0;
    virtual bool is_marked (const Node& node)  const = 0;

    virtual bool is_marked_branching (const Node& node) const = 0;

    virtual bool is_branching (const Node& node) const = 0;

    virtual void reset () = 0;

    virtual void dump () = 0;

protected:

    const Graph& _graph;
};


/********************************************************************************/
class NullTerminator :  public Terminator
{
public:

    static Terminator& singleton()
    {
        static Graph dummy;
        static NullTerminator instance(dummy); return instance;
    }

    virtual bool isEnabled () const { return false; }

    virtual void mark      (const Edge& edge) {}
    virtual bool is_marked (const Edge& edge)  const  { return false; };

    virtual void mark      (const Node& node) {}
    virtual bool is_marked (const Node& node)  const  { return false; }

    virtual bool is_marked_branching (const Node& node) const { return false; }

    virtual bool is_branching (const Node& node) const { return false; }

    virtual void reset () {}

    virtual void dump () {}

private:

    NullTerminator (const Graph& graph) : Terminator(graph)  {}
};

/********************************************************************************/
class BranchingTerminator :  public Terminator//, public ISmartIterator<Node>
{
public:

    typedef unsigned short int  Value;

    BranchingTerminator (const Graph& graph);
    ~BranchingTerminator();

    void mark      (const Edge& edge);
    bool is_marked (const Edge& edge)  const;

    void mark      (const Node& node);
    bool is_marked (const Node& node) const ;

    bool is_marked_branching (const Node& node) const ;

    bool is_branching (const Node& node) const ;

    void reset();

    void dump ();

private:

    bool is_indexed (const Node& node) const ;

    AssocSet<Node::Value, Value> branching_kmers;

    int getDelta (const Edge& edge) const;
};

/********************************************************************************/


#endif /* _GATB_TOOLS_TERMINATOR_HPP_ */

