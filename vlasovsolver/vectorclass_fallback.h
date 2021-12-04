/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#ifndef VECTORCLASS_PORTABLE_H
#define VECTORCLASS_PORTABLE_H
#include <math.h>
#include <iostream>
#include <vector>

/*! \file vectorclass_fallback.h
  \brief Simple class for implementing a vector with VECL real values

*/

// Prefetching does nothing in the fallback vectorclass, if no system implementation
// is available
#ifndef _mm_prefetch
#define _mm_prefetch(...)
#endif



template <class T>
class VecSimple {
public:
   T val[VECL] __attribute__((aligned(32)));
   // donot initi v
   VecSimple() { }
   // Replicate scalar x across v.

   template <typename ... S>
   VecSimple (S ... vargs){
      std::vector<T> vec = {vargs...};

      if (vec.size()==1){
         for (unsigned int i=0; i<VECL; i++){  
            val[i]=vec[0];
         }
      }else{
         for (unsigned int i=0; i<vec.size(); i++){  
            val[i]=vec[i];
         }
      }
   }





   // Copy vector v.
   VecSimple(VecSimple const &x){
      for(unsigned int i=0;i<VECL;i++)
         val[i]=x.val[i];
   }

   // Member function to load from array (unaligned)
   VecSimple & load(T const * p)  {
      for(unsigned int i=0;i<VECL;i++)
         val[i]=p[i];
      return *this;
   }
   // Member function to load from array, aligned by 32
   VecSimple & load_a(T const * p){
      return this->load(p);
   }
   
   VecSimple & insert(int i,T const &x) {
      val[i]=x;
      return *this;
   }


// Member function to store into array (unaligned)
   void store(T * p) const {
      for(unsigned int i=0;i<VECL;i++)
         p[i]=val[i];
   }
   // Member function to store into array, aligned by 32
   void store_a(T * p) const {
      this->store(p);
   }

   VecSimple & operator = (VecSimple const & r){
      for(unsigned int i=0;i<VECL;i++)
         val[i]=r.val[i];
      
      return *this;
   }

   T operator [](int i) const{
      return val[i];
   }

   VecSimple operator++ (int)
   {
      VecSimple<T> temp (*this);
      for(unsigned int i=0;i<VECL;i++)
         val[i]++;
      return temp;
   }
};

template <class T>
static inline VecSimple<T> abs(const VecSimple<T> &l){
   return VecSimple<T>(
      fabs(l.val[0]),
      fabs(l.val[1]),
      fabs(l.val[2]),
      fabs(l.val[3])
   );
}

template <class T>
static inline VecSimple<T> sqrt(const VecSimple<T> &l){
   return VecSimple<T>(
      sqrt(l.val[0]),
      sqrt(l.val[1]),
      sqrt(l.val[2]),
      sqrt(l.val[3])
   );
}



template <class T>
static inline VecSimple<T> operator + (const VecSimple<T> &l, const VecSimple<T> &r){
   return VecSimple<T>(
      l.val[0]+r.val[0],
      l.val[1]+r.val[1],
      l.val[2]+r.val[2],
      l.val[3]+r.val[3]
   );
}

template <class T, class S>
static inline VecSimple<T> operator + (const S &l, const VecSimple<T> &r){
   return VecSimple<T>(
      l+r.val[0],
      l+r.val[1],
      l+r.val[2],
      l+r.val[3]
   );
}

template <class T, class S>
static inline VecSimple<T> operator + (const VecSimple<T> &l, const S &r){
   return VecSimple<T>(
      l.val[0]+r,
      l.val[1]+r,
      l.val[2]+r,
      l.val[3]+r
   );
}
template <class T>
static inline VecSimple<T> operator - (const VecSimple<T> &r)
{
   return VecSimple<T>(
      -r.val[0],
      -r.val[1],
      -r.val[2],
      -r.val[3]
   );
}




template <class T>
static inline VecSimple<T> operator - (const VecSimple<T> &l, const VecSimple<T> &r)
{
   return VecSimple<T>(
      l.val[0]-r.val[0],
      l.val[1]-r.val[1],
      l.val[2]-r.val[2],
      l.val[3]-r.val[3]
   );
}

template <class T, class S>
static inline VecSimple<T> operator - (const S &l, const VecSimple<T> &r){
   return VecSimple<T>(
      l-r.val[0],
      l-r.val[1],
      l-r.val[2],
      l-r.val[3]
   );
}

template <class T, class S>
static inline VecSimple<T> operator - (const VecSimple<T> &l, const S &r){
   return VecSimple<T>(
      l.val[0]-r,
      l.val[1]-r,
      l.val[2]-r,
      l.val[3]-r
   );
}

template <class T>
static inline VecSimple<T> operator * (const VecSimple<T> &l, const VecSimple<T> &r)
{
   return VecSimple<T>(
      l.val[0]*r.val[0],
      l.val[1]*r.val[1],
      l.val[2]*r.val[2],
      l.val[3]*r.val[3]
   );
}

template <class T, class S>
static inline VecSimple<T> operator * (const VecSimple<T> &l, const S &r)
{
   return VecSimple<T>(
      l.val[0]*r,
      l.val[1]*r,
      l.val[2]*r,
      l.val[3]*r
   );
}

template <class T, class S>
static inline VecSimple<T> operator * (const S &l,const VecSimple<T> &r)
{
   return VecSimple<T>(
      l*r.val[0],
      l*r.val[1],
      l*r.val[2],
      l*r.val[3]
   );
}



template <class T>
static inline VecSimple<T> operator / (const VecSimple<T> &l, const VecSimple<T> &r)
{
   return VecSimple<T>(
      l.val[0]/r.val[0],
      l.val[1]/r.val[1],
      l.val[2]/r.val[2],
      l.val[3]/r.val[3]
   );
}

template <class T, class S>
static inline VecSimple<T> operator / (const VecSimple<T> &l, const S &r)
{
   return VecSimple<T>(
      l.val[0]/r,
      l.val[1]/r,
      l.val[2]/r,
      l.val[3]/r
   );
}

template <class T, class S>
static inline VecSimple<T> operator / (const S &l, const VecSimple<T> &r )
{
   return VecSimple<T>(
      l/r.val[0],
      l/r.val[1],
      l/r.val[2],
      l/r.val[3]
   );
}

template <class T>
static inline  VecSimple<T> & operator += (VecSimple<T> &l, const VecSimple<T> &r){
   l=l+r;
   return l;
}

template <class T, class S>
static inline  VecSimple<T> & operator += (VecSimple<T> &l, const S &r){
   l = l+r;
   return l;
}

template <class T>
static inline VecSimple<T> & operator -= (VecSimple<T> &l, const VecSimple<T> &r){
   l=l-r;
   return l;
}

template <class T, class S>   
static inline VecSimple<T> & operator -= (VecSimple<T> &l, const S &r){
   l = l - r;
   return l;
}

template <class T>
static inline VecSimple<bool> operator || (const VecSimple<T> &l, const VecSimple<T> &r)
{
   return VecSimple<bool>(
      l.val[0] || r.val[0],
      l.val[1] || r.val[1],
      l.val[2] || r.val[2],
      l.val[3] || r.val[3]
   );
}



template <class T>
static inline VecSimple<bool> operator && (const VecSimple<T> &l, const VecSimple<T> &r)
{
   return VecSimple<bool>(
      l.val[0] && r.val[0],
      l.val[1] && r.val[1],
      l.val[2] && r.val[2],
      l.val[3] && r.val[3]
   );
}

template <class T>
static inline VecSimple<bool> operator == (const VecSimple<T> &l, const VecSimple<T> &r)
{
   return VecSimple<bool>(
      l.val[0] == r.val[0],
      l.val[1] == r.val[1],
      l.val[2] == r.val[2],
      l.val[3] == r.val[3]
   );
}

template <class T, class S>
static inline VecSimple<bool> operator == (const VecSimple<T> &l, const S& r)
{
   return VecSimple<bool>(
      l.val[0] == r,
      l.val[1] == r,
      l.val[2] == r,
      l.val[3] == r
   );
}

template <class T, class S>
static inline VecSimple<bool> operator != (const VecSimple<T> &l, const S& r)
{
   return VecSimple<bool>(
      l.val[0] != r,
      l.val[1] != r,
      l.val[2] != r,
      l.val[3] != r
   );
}

template <class T>
static inline VecSimple<bool> operator ! (const VecSimple<T> &l)
{
   return VecSimple<bool>(
      !l.val[0],
      !l.val[1],
      !l.val[2],
      !l.val[3]
   );
}


template <class T>
static inline VecSimple<bool> operator > (const VecSimple<T> &l, const VecSimple<T> &r)
{
   return VecSimple<bool>(
      l.val[0] > r.val[0],
      l.val[1] > r.val[1],
      l.val[2] > r.val[2],
      l.val[3] > r.val[3]
   );
}


template <class T, class S>
static inline VecSimple<bool> operator > (const VecSimple<T> &l, const S r)
{
   return VecSimple<bool>(
      l.val[0] > r,
      l.val[1] > r,
      l.val[2] > r,
      l.val[3] > r
   );
}



template <class T, class S>
  static inline VecSimple<bool> operator > (const S l,const VecSimple<T> &r) 
{
   return VecSimple<bool>(
      l > r.val[0],
      l > r.val[1],
      l > r.val[2],
      l > r.val[3]
   );
}


template <class T>
static inline VecSimple<bool> operator >= (const VecSimple<T> &l, const VecSimple<T> &r)
{
   return VecSimple<bool>(
      l.val[0] >= r.val[0],
      l.val[1] >= r.val[1],
      l.val[2] >= r.val[2],
      l.val[3] >= r.val[3]
   );
}


template <class T, class S>
static inline VecSimple<bool> operator >= (const VecSimple<T> &l, const S r)
{
   return VecSimple<bool>(
      l.val[0] >= r,
      l.val[1] >= r,
      l.val[2] >= r,
      l.val[3] >= r
   );
}



template <class T, class S>
  static inline VecSimple<bool> operator >= (const S l,const VecSimple<T> &r) 
{
   return VecSimple<bool>(
      l >= r.val[0],
      l >= r.val[1],
      l >= r.val[2],
      l >= r.val[3]
   );
}



template <class T>
static inline VecSimple<bool> operator < (const VecSimple<T> &l, const VecSimple<T> &r)
{
   return VecSimple<bool>(
      l.val[0] < r.val[0],
      l.val[1] < r.val[1],
      l.val[2] < r.val[2],
      l.val[3] < r.val[3]
   );
}


template <class T, class S>
static inline VecSimple<bool> operator < (const VecSimple<T> &l,const S &r)
{
   return VecSimple<bool>(
      l.val[0] < r,
      l.val[1] < r,
      l.val[2] < r,
      l.val[3] < r
   );
}

template <class T, class S>
  static inline VecSimple<bool> operator < (const S l, const VecSimple<T> &r) 
{
   return VecSimple<bool>(
      l < r.val[0],
      l < r.val[1],
      l < r.val[2],
      l < r.val[3]
   );
}



template <class T>
static inline VecSimple<bool> operator <= (const VecSimple<T> &l, const VecSimple<T> &r)
{
   return VecSimple<bool>(
      l.val[0] <= r.val[0],
      l.val[1] <= r.val[1],
      l.val[2] <= r.val[2],
      l.val[3] <= r.val[3]
   );
}


template <class T, class S>
static inline VecSimple<bool> operator <= (const VecSimple<T> &l,const S &r)
{
   return VecSimple<bool>(
      l.val[0] <= r,
      l.val[1] <= r,
      l.val[2] <= r,
      l.val[3] <= r
   );
}

template <class T, class S>
  static inline VecSimple<bool> operator <= (const S l, const VecSimple<T> &r) 
{
   return VecSimple<bool>(
      l <= r.val[0],
      l <= r.val[1],
      l <= r.val[2],
      l <= r.val[3]
   );
}




template <class T>
static inline VecSimple<T> min(VecSimple<T> const & l, VecSimple<T> const & r){
   return VecSimple<T>(
      l.val[0]<r.val[0]?l.val[0]:r.val[0],
      l.val[1]<r.val[1]?l.val[1]:r.val[1],
      l.val[2]<r.val[2]?l.val[2]:r.val[2],
      l.val[3]<r.val[3]?l.val[3]:r.val[3]
   );
}

template <class T, class S>
static inline VecSimple<T> min(S const & l, VecSimple<T> const & r){
   return VecSimple<T>(
      l <r.val[0] ? l:r.val[0],
      l <r.val[1] ? l:r.val[1],
      l <r.val[2] ? l:r.val[2],
      l <r.val[3] ? l:r.val[3]
   );
}

template <class T>
static inline VecSimple<T> max(VecSimple<T> const & l, VecSimple<T> const & r){
   return VecSimple<T>(
      l.val[0]>r.val[0]?l.val[0]:r.val[0],
      l.val[1]>r.val[1]?l.val[1]:r.val[1],
      l.val[2]>r.val[2]?l.val[2]:r.val[2],
      l.val[3]>r.val[3]?l.val[3]:r.val[3]
   );
}


template <class T, class S>
static inline VecSimple<T> max(VecSimple<T> const & l, S const & r){
   return VecSimple<T>(
      l.val[0] > r ? l.val[0] : r,
      l.val[1] > r ? l.val[1] : r,
      l.val[2] > r ? l.val[2] : r,
      l.val[3] > r ? l.val[3] : r
   );
}


template <class T, class S>
  static inline VecSimple<T> max(S const & l, VecSimple<T> const & r){
  return VecSimple<T>(
     r.val[0] > l ? r.val[0] : l,
     r.val[1] > l ? r.val[1] : l,
     r.val[2] > l ? r.val[2] : l,
     r.val[3] > l ? r.val[3] : l
  );
}



template <class T>
static inline VecSimple<T> select(VecSimple<bool> const & a, VecSimple<T> const & b, VecSimple<T> const & c){ 
   return VecSimple<T>(
      a.val[0] ? b.val[0] : c.val[0],
      a.val[1] ? b.val[1] : c.val[1],
      a.val[2] ? b.val[2] : c.val[2],
      a.val[3] ? b.val[3] : c.val[3]
   );
}


template <class T, class S>
static inline VecSimple<T> select(VecSimple<bool> const & a, S const & b, VecSimple<T> const & c){ 
   return VecSimple<T>(
      a.val[0] ? b : c.val[0],
      a.val[1] ? b : c.val[1],
      a.val[2] ? b : c.val[2],
      a.val[3] ? b : c.val[3]
   );
}


template <class T, class S>
  static inline VecSimple<T> select(VecSimple<bool> const & a, VecSimple<T> const & b, S const & c){
   return VecSimple<T>(
      a.val[0] ? b.val[0] : c,
      a.val[1] ? b.val[1] : c,
      a.val[2] ? b.val[2] : c,
      a.val[3] ? b.val[3] : c
   );
}


template <class T>
  static inline VecSimple<T> select(VecSimple<bool> const & a, T const & b, T const & c){
   return VecSimple<T>(
      a.val[0] ? b : c,
      a.val[1] ? b : c,
      a.val[2] ? b : c,
      a.val[3] ? b : c
   );
}

template <class T>
static inline bool horizontal_or(VecSimple<T> const & a){ 
  return a.val[0] || a.val[1] || a.val[2] || a.val[3];
}


template <class T>
static inline bool horizontal_and(VecSimple<T> const & a){ 
  return a.val[0] && a.val[1] && a.val[2] && a.val[3];
}



template <class T>
static inline VecSimple<int> truncate_to_int(VecSimple<T> const & a){ 
  return VecSimple<int>(a.val[0], a.val[1], a.val[2], a.val[3]);
}

template <class T>
static inline VecSimple<double> to_double(VecSimple<T> const & a){ 
  return VecSimple<double>(a.val[0], a.val[1], a.val[2], a.val[3]);
}

template <class T>
static inline VecSimple<float> to_float(VecSimple<T> const & a){ 
  return VecSimple<float>(a.val[0], a.val[1], a.val[2], a.val[3]);
}

// Dummy functions that Agner vectorclass has.
// These are here to suppress compiler error messages only
static void no_subnormals() 
{
}





#endif
