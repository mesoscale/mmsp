#ifndef priority_queue_h_
#define priority_queue_h_

#include <iostream>
#include <vector>
#include <map>
#include <cassert>

// The DistanceVoxel_PriorityQueue is a customized, non-templated
// priority queue that stores DistanceVoxel pointers in a heap.  The
// elements in the heap can be looked up in a map, to quickly find out
// the current index of the element within the heap.

// =========================================================================

class DistanceVoxel_PriorityQueue
{

public:
  // --------------------------
  // CONSTRUCTORS
  DistanceVoxel_PriorityQueue(){} // Default Constructor
  DistanceVoxel_PriorityQueue( const std::vector<DistanceVoxel* > &values )
  { // Class Constructor
    // Build heap from input vector using PUSH member function
    for ( std::vector<DistanceVoxel* >::const_iterator itr = values.begin(); itr != values.end(); ++itr)
    {
      push( *itr );
    }
  }
  // NOTE: No Copy Constructor or Assignment Operator are defined.

  // ------------------------
  // ACCESSORS
  int size()
  {
    return m_heap.size();
  }
  bool empty()
  {
    return m_heap.empty();
  }
  int last_non_leaf()
  {
    return ( size() - 1 ) / 2;
  }
  int get_parent( int i )
  {
    assert ( i > 0 );
    assert( i < size() );
    return ( i - 1 ) / 2;
  }
  int has_left_child( int i )
  {
    return ( 2 * i ) + 1 < size();
  }
  int has_right_child( int i )
  {
    return ( 2 * i ) + 2 < size();
  }
  int get_left_child( int i )
  {
    assert ( i >= 0 && has_left_child( i ) );
    return 2 * i + 1;
  }
  int get_right_child( int i )
  {
    assert ( i >= 0 && has_right_child( i ) );
    return 2 * i + 2;
  }

  // read the top element
  const DistanceVoxel* top() const
  {
    assert( !m_heap.empty() );
    return m_heap[0];
  }

  // is this element in the heap?
  bool in_heap( DistanceVoxel* element ) const
  {
    std::map<DistanceVoxel*, int>::const_iterator itr = backpointers.find( element );
    return ( itr != backpointers.end() );
  }

  // add an element to the heap
  void push( DistanceVoxel* element )
  {
    std::map<DistanceVoxel*, int>::iterator itr = backpointers.find( element );
    assert ( itr == backpointers.end() );
    m_heap.push_back( element );
    backpointers[element] = m_heap.size() - 1;
    this->percolate_up( int( m_heap.size() - 1 ) );
  }

  // the value of this element has been edited, move the element up or down
  void update_position( DistanceVoxel* element )
  {
    std::map<DistanceVoxel*, int>::iterator itr = backpointers.find( element );
    assert ( itr != backpointers.end() );
    this->percolate_up( itr->second );
    this->percolate_down( itr->second );
  }

  // remove the top (minimum) element
  void pop()
  {
    assert( !m_heap.empty() );
    int success = backpointers.erase( m_heap[0] );
    assert ( success == 1 );
    m_heap[0] = m_heap.back();
    m_heap.pop_back();
    this->percolate_down( 0 );
  }

private:
  // REPRESENTATION
  //  the heap is stored in a vector representation (the binary tree
  //  structure "unrolled" one row at a time)
  std::vector<DistanceVoxel*> m_heap;

  // the map stores a correpondence between elements & indices in the heap
  std::map<DistanceVoxel*, int> backpointers;

  // Helper function for percolators: swap contents of heap
  void swap( int i, int j)
  {
    // Change positions in heap
    DistanceVoxel* temp = m_heap[j];
    m_heap[j] = m_heap[i];
    m_heap[i] = temp;
    // Update positions in map
    backpointers[m_heap[i]] = i;
    backpointers[m_heap[j]] = j;
  }

  // private helper functions
  void percolate_up( int i )
  {
    while ( i > 0 )
    { // Loop until you become the root node
      if ( m_heap[i]->getValue() < m_heap[get_parent(i)]->getValue() )
      {
        swap(i, get_parent(i));
        i = get_parent(i);
      }
      else
      {
        break;
      }
    }
  }

  void percolate_down( int i )
  {
    while ( has_left_child(i) )
    { // Loop as long as there's at least one level below you
      int child;
      // Choose the child to compare against
      if ( has_right_child(i) && m_heap[ get_right_child(i) ]->getValue() < m_heap[ get_left_child(i) ]->getValue() )
      {
        child = get_right_child(i);
      }
      else child = get_left_child(i);
      if ( m_heap[child]->getValue() < m_heap[i]->getValue() )
      {
        swap(child,i);
        i = child;
      }
      else
      {
        break;
      }
    }
  }

};

#endif
