/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#ifndef RANGE_TREE_HPP
#define RANGE_TREE_HPP

#include <set>

template <class T>
class RangeTree
{
  public:
  
    struct Span 
    {
      mutable T first, last;
      bool operator<( const Span& oth ) const { return last < oth.first; }
    };
 
  private:
  
    typedef std::set<Span> SpanTree;
    SpanTree set_;
  
  public:
  
    class const_iterator
    {
      public:
      
        const_iterator()             : value_(0) {}
        
        const_iterator( typename SpanTree::const_iterator span, 
                        typename SpanTree::const_iterator end,
                        T value )    : beg_(span), span_(span), end_(end), value_(value) {}
        
        T operator*() const          { return value_; }
        
        const_iterator& operator++();
        const_iterator& operator--();
        
        const_iterator  operator++(int)
          { const_iterator tmp(*this); operator++(); return tmp; }
          
        const_iterator  operator--(int)
          { const_iterator tmp(*this); operator--(); return tmp; }
        
        bool operator==( const const_iterator& other ) const
          { return span_ == other.span_ && value_ == other.value_; }
        
        bool operator!=( const const_iterator& other ) const
          { return span_ != other.span_ || value_ != other.value_; }
        
      
      private:
      
        typename SpanTree::const_iterator beg_;
        typename SpanTree::const_iterator span_;
        typename SpanTree::const_iterator end_;
        T value_;
    };
    
    typedef const_iterator iterator;
    
    const_iterator begin() const 
      { return const_iterator(set_.begin(), set_.end(), set_.begin()->first); }
    
    const_iterator end() const
      { return const_iterator(set_.end(), set_.end(), 0); }
    
    const_iterator find( T value ) const;
    
    const_iterator insert( T value );
    
    T size() const;
    
    long num_ranges() const { return set_.size(); }
    
    bool empty() const;
    
    class span_iterator
    {
      private:
        typename SpanTree::const_iterator iter_;
      
      public:
      
        span_iterator( ) {}
        span_iterator( typename SpanTree::const_iterator iter ) 
          : iter_(iter) {}
        
        std::pair<T,T> operator*() const
          { const std::pair<T,T> tmp( iter_->first, iter_->last ); return tmp; }

        span_iterator& operator++() 
          { iter_++; return *this; }
          
        span_iterator& operator--()
          { iter_--; return *this; }
        
        span_iterator  operator++(int)
          { return span_iterator( iter_++ ); }
          
        span_iterator  operator--(int)
          { return span_iterator( iter_-- ); }
          
        bool operator==( const span_iterator& other )
          { return iter_ == other.iter_; }
          
        bool operator!=( const span_iterator& other )
          { return iter_!= other.iter_; }
    };
    
    span_iterator span_begin() const
      { return span_iterator( set_.begin() ); }
    
    span_iterator span_end() const
      { return span_iterator( set_.end() ); }
};

#include "RangeTree.cpp"

#endif
