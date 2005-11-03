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
  
    typedef std::set<Span> SpanTree;
 
  private:
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


template <class T> typename RangeTree<T>::const_iterator& 
RangeTree<T>::const_iterator::operator++()
{
  if (span_ == end_)
    return *this;
  
  if (value_ == span_->last)
  {
    ++span_;
    value_ = (span_ == end_) ? 0 : span_->first;
  }
  else
  { 
    ++value_;
  }
  
  return *this;
}


template <class T> typename RangeTree<T>::const_iterator& 
RangeTree<T>::const_iterator::operator--()
{
  if (span_ == end_)
  {
    --span_;
    value_ = span_->last;
  }
  
  else if (value_ == span_->first && span_ != beg_)
  {
    --span_;
    value_ = span_->first;
  }
  
  else
  { 
    --value_;
  }
  
  return *this;
}

template <class T> typename RangeTree<T>::const_iterator
RangeTree<T>::find( T value ) const
{
  Span span = { value, value };
  typename SpanTree::iterator iter = set_.find( span );
  if (iter == set_.end())
    return end();
  else
    return const_iterator( iter, set_.end(), value );  
}

template <class T> typename RangeTree<T>::const_iterator
RangeTree<T>::insert( T value )
{
  Span span = { value, value };
  typename SpanTree::iterator iter = set_.find( span );
  if (iter != set_.end()) 
    return const_iterator(iter, set_.end(), value);
  
  --span.first;
  iter = set_.find( span );
  if (iter != set_.end())
  {
    typename SpanTree::iterator next = iter;
    ++next;
    if (next != set_.end() && next->first - 1 == value )
    {
      T tmp = next->last;
      set_.erase( next );
      iter->last = tmp;
    }
    else
    {
      ++iter->last;
    }
    return const_iterator(iter, set_.end(), value);
  }
  ++span.first;
  
  ++span.last;
  iter = set_.find( span );
  if (iter != set_.end())
  {
    iter->first--;
    return const_iterator(iter, set_.end(), value);
  }
  --span.last;
  
  iter = set_.insert( span ).first;
  return const_iterator(iter, set_.end(), value);
}

template <class T> T RangeTree<T>::size( ) const
{
  T count = 0;
  for (typename SpanTree::iterator iter = set_.begin(); iter != set_.end(); ++iter)
    count += iter->last - iter->first + 1;
  return count;
}

template <class T> bool RangeTree<T>::empty() const
{
  return set_.empty();
}


#endif
