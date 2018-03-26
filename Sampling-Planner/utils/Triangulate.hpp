

#ifndef M_TRIANGULATE
#define M_TRIANGULATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>  // Include STL vector class.

class m_Vector2d
{
public:
  m_Vector2d(float x,float y)
  {
    Set(x,y);
  };

  float GetX(void) const { return mX; };

  float GetY(void) const { return mY; };

  void  Set(float x,float y)
  {
    mX = x;
    mY = y;
  };
private:
  float mX;
  float mY;
};

// Typedef an STL vector of vertices which are used to represent
// a polygon/contour and a series of triangles.
typedef std::vector< m_Vector2d > Vector2dVector;


class m_Triangulate
{
public:

  // triangulate a contour/polygon, places results in STL vector
  // as series of triangles.
  static bool Process(const Vector2dVector &contour,
                      Vector2dVector &result);

  // compute area of a contour/polygon
  static float Area(const Vector2dVector &contour);

  // decide if point Px/Py is inside triangle defined by
  // (Ax,Ay) (Bx,By) (Cx,Cy)
  static bool InsideTriangle(float Ax, float Ay,
                      float Bx, float By,
                      float Cx, float Cy,
                      float Px, float Py);


private:
  static bool Snip(const Vector2dVector &contour,int u,int v,int w,int n,int *V);

};

#endif // M_TRIANGULATE

