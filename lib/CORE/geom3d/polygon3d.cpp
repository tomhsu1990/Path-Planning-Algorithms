/*****************************************************************
 * File: polygon3d.cc
 * Synopsis:
 *      Basic 3-dimensional geometry
 * Author: Shubin Zhao (shubinz@cs.nyu.edu), 2001.
 *
 *****************************************************************
 * CORE Library Version 1.4 (July 2001)
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * Copyright (c) 1995, 1996, 1998, 1999, 2000, 2001 Exact Computation Project
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id: polygon3d.cpp,v 1.1 2006/04/03 18:55:31 exact Exp $
 *****************************************************************/

#include <CORE/geom3d/polygon3d.h>

/************************************************************
  *   constructors
************************************************************/

Polygon3d::Polygon3d() {
  headN = NULL;
  size = 0;
}

  // copy constructor
Polygon3d::Polygon3d(const Polygon3d& plg) {
  copy( plg );
  size = plg.getSize();
}

Polygon3d::~Polygon3d() {
  freeMemory(); 
}

/************************************************************
  *   member functions
************************************************************/
void Polygon3d::freeMemory() {
  if( headN != NULL ) // free memory
    for( int i=0; i<size; i++ ) {
      PointNode* pNode = headN;
      headN = headN->next;
      delete pNode;
    }

  headN = NULL;
  size = 0;
}

  // copy point list from other polygon
  // make sure old list has been deleted before calling this
void Polygon3d::copy( const Polygon3d& plg) {
  headN = NULL;
  size = 0;
  for( int i=0; i<plg.getSize(); i++ ) 
    addPoint( *plg[i] );
}

bool Polygon3d::searchPoint( const Point3d& p ) const {
  PointNode* pNode = headN;
  for( int i=0; i<size; i++, pNode = pNode->next )
    if( p == *(pNode->p) )
      return true;

  return false;
}

Point3d* Polygon3d::getPoint( int index ) const {
  assert( index < size );
  PointNode* pNode = headN;
   //locate
  for( int i=0; i<index; i++, pNode = pNode->next );
  return pNode->p;
}

void Polygon3d::removePoint( int index ) {
  assert( index < size );
  PointNode* pNode = headN;
  if( index == 0 )  headN = headN->next;

   //locate
  for( int i=0; i<index; i++, pNode = pNode->next );
  pNode->prev->next = pNode->next;
  pNode->next->prev = pNode->prev;
  delete pNode;
  size--;

  if( size==0 ) headN = NULL;
}

bool Polygon3d::removePoint( const Point3d& p ){
  PointNode* pNode = headN;
  if( *(pNode->p) == p ) // move list head 
    headN = pNode->next;

    for( int i=0; i<size; i++ ) {
      if( *(pNode->p) == p )  {
        pNode->prev->next = pNode->next;
        pNode->next->prev = pNode->prev;
        delete pNode;
        size--;
        if( size == 0 )  headN = NULL;   
        return true;
      }

    pNode = pNode->next;
  }

  return false;
}

void Polygon3d::removeAllPoints() {
  freeMemory();
}

 // add a new point
void Polygon3d::addPoint( const Point3d& p  ) {
  if( searchPoint(p) )  return;  // do not add duplicated point

  if( headN == NULL ) {
    headN = new PointNode(p);
    headN->next = headN;
    headN->prev = headN;
  }
  else {  //insert in a double linked list
    PointNode *pNode = new PointNode(p);
    headN->prev->next = pNode;
    pNode->next = headN;
    pNode->prev = headN->prev;
    headN->prev = pNode;
  }

  size++;
}

 // get previous point of given point
Point3d* Polygon3d::nextPoint( const Point3d& pt )  const { 
  PointNode* pNode = headN;
  for( int i=0; i<size; i++, pNode=pNode->next ) 
    if( pt == *(pNode->p) ) 
      return pNode->next->p;

   // p is not found
  return NULL;
}

 // get previous point of given point
Point3d* Polygon3d::prevPoint( const Point3d& pt ) const  { 
  PointNode* pNode = headN;
  for( int i=0; i<size; i++, pNode=pNode->next ) 
    if( pt == *(pNode->p) ) 
      return pNode->prev->p;

   // p is not found
  return NULL;
}

Polygon3d& Polygon3d::operator=( const Polygon3d& plg)
{
  removeAllPoints();
  copy( plg );
  size = plg.getSize();
  return *this;
} 

 // test if two polygon is identical
 // two polygons are indentical iff they have the same set of points 
 // in the same order
bool Polygon3d::operator ==(Polygon3d& plg) const {
  if( size != plg.getSize() )  return false;
  assert( size >= 2);

   // locate the first point of plg
  Iterator I = plg.getIterator();
  Point3d * pPoint = I.getPoint();
  PointNode * pNode = headN;
  int i;
  for( i=0; i<size; i++, pNode=pNode->next ) 
    if( *(pNode->p) == *pPoint ) break;

  if( i == size )   // not found the first point
    return false;

   // we should check both directions
  Iterator J = I;
  PointNode * pTempN = pNode;
  for( i=0; i<size; i++, J++, pTempN=pTempN->next ) 
    if( *(J.getPoint()) != *(pTempN->p) ) break;
 
  if( i==size ) return true;

  for( i=0, J=I, pTempN = pNode; i<size; i++, J--, pTempN=pTempN->next ) 
    if( *(J.getPoint()) != *(pTempN->p) ) break;
 
  if( i==size ) return true;

  return false;

}

bool Polygon3d::verify() {
  if( size <= 3 ) {
    std::cout<< "\nWARNING:Polygon has less than three points";
    return false;
  }

  // 1. verify coplanarity
  // get container plane first
  Plane3d pl = Plane3d( *getPoint(0), *getPoint(1), *getPoint(2) ); 

  PointNode* pNode = headN;
  for( int i=0; i<size; i++ ) {
    if( !pl.contains( *(pNode->p) ) ) 
      return false;

    pNode = pNode->next;
  }

  // 2. verify edges
  // to be compeled

  return true;
}
   
   // test coplanarity 
bool Polygon3d::isCoplanar( const Point3d& p ) const {
  assert( size>= 3 );
  Plane3d pl = Plane3d( *getPoint(0), *getPoint(1), *getPoint(2) );
  return pl.contains( p );
}  

   // test coplanarity 
bool Polygon3d::isCoplanar( const Plane3d& pl ) const {
  PointNode* pNode = headN;
  for( int i=0; i<size; i++ ) {
    if( !pl.contains( *(pNode->p) ) ) 
      return false;

    pNode = pNode->next;
  }

  return true;
}


   // test if p is on amy edge
bool Polygon3d::isOnEdge( const Point3d& p ) const {
  PointNode* pNode = headN;
  for( int i=0; i<size; i++ ) {
    Segment3d s( *(pNode->p), *(pNode->next->p) );
    if( s.contains(p) ) 
      return true;

    pNode = pNode->next;
  }

  return false;
}


/************************************************************
  *   I/O 
 ************************************************************/

std::ostream &operator<<(std::ostream &o, const Polygon3d &plg) {
  o << "Polygon3d: size=" << plg.getSize() << std::endl;
  for( int i=0; i<plg.getSize(); i++ )
    std::cout << "\t" << *plg[i] << std::endl;
  return o;
}
