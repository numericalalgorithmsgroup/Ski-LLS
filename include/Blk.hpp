/**
 * @file   blk.hpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Thu Jul 28 15:10:23 2011
 * 
 * @brief  memory block
 * 
 * 
 */
#ifndef _BLK_HPP
#define _BLK_HPP

#include <cstddef>
#include "boost/shared_array.hpp"

#include "Config.hpp"

template <typename T>
class Blk;

typedef Blk<double> Blk_d;
typedef Blk<long>   Blk_l;

/**
 * A memory block class using boost::shared_array to deallocate automatically.
 * 
 */ 
template <typename T>
class Blk
{
 public:

  typedef T      val_t;   
  typedef Blk<T> self_t;  

  /** 
   * constructor
   * 
   * @param size blk size
   */
  Blk( std::size_t size = 0 )  
      : _sz(size), _ptr(new T[size]) 
  {}

  T* begin() const { return _ptr.get(); } 
  T* end() const { return _ptr.get() + _sz; } 

  std::size_t size() const { return _sz; } 
  bool empty() const { return _sz==0; } 

  bool has( const T* ptr ) const
  {
    return ( ptr >= begin() ) && ( ptr < end() ); 
  }

  const T operator()( std::size_t i ) const { return *(_ptr.get()+i); } 
  T& operator()( std::size_t i ) { return *(_ptr.get()+i); } 
  
 private:
  
  std::size_t _sz; 
  boost::shared_array<T> _ptr; 
};

#endif  // _BLK_HPP
