#include <CORE/geom2d/polygon2d.h>

Polygon2d::Polygon2d(const std::vector<Point2d>& Pts){
  for(unsigned i=0;i<Pts.size();++i){
    pts.push_back(Pts[i]);
  }
  cp = center();
  orient = clockwise();
}

Polygon2d::~Polygon2d(){
}

Polygon2d::Polygon2d(){}

Polygon2d::Polygon2d(const Polygon2d& poly){
  for(unsigned i=0;i<poly.pts.size();++i){
    pts.push_back(poly.pts[i]);
  }
  cp = poly.cp;
  orient = poly.orient;
}

Polygon2d& Polygon2d::operator=(const Polygon2d& poly){
  for(unsigned i=0;i<poly.pts.size();++i){
    pts.push_back(poly.pts[i]);
  }
  cp = poly.cp;
  orient = poly.orient;
  return *this;
} 

Point2d Polygon2d::center() {
  cp = Point2d(0,0);
  double cnt = 0;
  for(unsigned i=0;i<pts.size();++i){
    cp.set(cp.X()+pts[i].X(), cp.Y()+pts[i].Y());
    ++cnt;
  }
  if(cnt > 0) cp.set(cp.X()/cnt, cp.Y()/cnt);

  return cp;
}

bool Polygon2d::operator==(const Polygon2d& poly) const {
  for(unsigned i=0;i<pts.size();++i){
    if(!poly.contains(pts[i])) return false;
  }
  return true;
}

double Polygon2d::distance(const Point2d& p){
  double minDist = pts[0].distance(p);
  for(unsigned i=1;i<pts.size()-1;++i){
    double dist = pts[i].distance(p);
    if(minDist > dist){
      minDist = dist;
    }
  }
  for(unsigned i=0;i<pts.size()-1;++i){
    Segment2d seg(pts[i], pts[i+1]);
    double dist = seg.distance(p);
    if(minDist > dist){
      minDist = dist;
    }
  }
  return minDist;
}

double Polygon2d::distance(const Line2d& l){
  double minDist = l.distance(pts[0]);
  for(unsigned i=1;i<pts.size()-1;++i){
    double dist = l.distance(pts[i]);
    if(minDist > dist){
      minDist = dist;
    }
  }
  return minDist;
}

double Polygon2d::distance(Polygon2d& poly) {
   return 0;
}

int Polygon2d::side_of(const Point2d& p) const {
  for(unsigned i=0;i<pts.size()-1;++i){
    Segment2d seg(pts[i], pts[i+1]);
    if(seg.contains(p)){
      return 0;
    }
  }
  // assume the orientation of the polygon is clockwise
  for(unsigned i=0;i<pts.size()-1;++i){
    Segment2d seg(pts[i], pts[i+1]);
    if(seg.orientation(p) < 0){
      return -1;
    }
  }
  return 1;
} 
 
bool Polygon2d::outside(const Point2d& p) {
   return (orient * side_of(p)) < 0;
}

bool Polygon2d::inside(const Point2d& p) {

    if (pts.size() < 3)  return false;

    // Create a point for line segment from p to infinite
    Point2d extreme(INT_MAX, p.Y());
    Segment2d seg_p_ext(p, extreme);

    // Count intersections of the above line with sides of polygon
    int count = 0, i = 0;
    do {
        int next = (i+1)%pts.size();
        Segment2d seg(pts[i], pts[next]);
        // Check if the line segment from 'p' to 'extreme' intersects
        // with the line segment seg
        if (seg.intersects(seg_p_ext) >= 0) {
            // If the point 'p' is colinear with line segment 'i-next',
            // then check if it lies on segment. If it lies, return true,
            // otherwise false
            if (seg.onSegment(p)) return true;

            ++count;
        }
        i = next;
    } while (i != 0);

    // Return true if count is odd, false otherwise
    return count&1;  // Same as (count%2 == 1)

   // only for convex polygon
   //return (orient * side_of(p)) >= 0;
}

bool Polygon2d::contains(const Point2d& p) const {
   return side_of(p) == 0;
}

std::ostream& operator<<(std::ostream& out, Polygon2d& poly) {
   out << "Polygon2d[ center=" << poly.center() << "]";
   return out;
}

std::istream& operator>>(std::istream& in, Polygon2d poly) {
  std::vector<Point2d> pts;
  Point2d tmp;
  int i=1;
  do{
      std::cout << "\nInput " << i << " point(-,-):";
      if(!(in >> tmp)) break;
      ++i;
      pts.push_back(tmp);
  }while(true);
  poly = Polygon2d(pts);
  return in;
}

int Polygon2d::clockwise(){
  if(pts.size() == 0){
    return 0;
  }

  for(unsigned i=0;i<pts.size()-2;++i){
    Segment2d seg(pts[i], pts[i+1]);
    if(seg.orientation(pts[i+2]) < 0){
      return -1;
    }
  }
  return 1;
}

float Polygon2d::area(const std::vector<Point2d> &contour){

  int n = contour.size();
  float A=0.0f;

  for(int p=n-1,q=0; q<n; p=q++){
    A += contour[p].X()*contour[q].Y() - contour[q].X()*contour[p].Y();
  }
  return A*0.5f;
}

bool Polygon2d::processTriangulate(const std::vector<Point2d> &contour, std::vector<Point2d> &result){
  /* allocate and initialize list of Vertices in polygon */

  int n=contour.size();
  if(n<3) return false;

  int *V = new int[n];

  /* we want a counter-clockwise polygon in V */
  if(0.0f < area(contour))
    for(int v=0; v<n; ++v) V[v] = v;
  else
    for(int v=0; v<n; ++v) V[v] = (n-1)-v;

  int nv = n;

  /*  remove nv-2 Vertices, creating 1 triangle every time */
  int count = 2*nv;   /* error detection */

  for(int m=0, v=nv-1; nv>2;){
    /* if we loop, it is probably a non-simple polygon */
    if(0 >= (count--)){
      //** Triangulate: ERROR - probable bad polygon!
      return false;
    }

    /* three consecutive vertices in current polygon, <u,v,w> */
    int u = v  ; if (nv <= u) u = 0; /* previous */
    v = u+1;     if (nv <= v) v = 0; /* new v    */
    int w = v+1; if (nv <= w) w = 0; /* next     */

    if(snip(contour,u,v,w,nv,V)){
      int a,b,c,s,t;

      /* true names of the vertices */
      a = V[u]; b = V[v]; c = V[w];

      /* output Triangle */
      result.push_back( contour[a] );
      result.push_back( contour[b] );
      result.push_back( contour[c] );

      m++;

      /* remove v from remaining polygon */
      for(s=v,t=v+1;t<nv;s++,t++) V[s] = V[t]; nv--;

      /* resest error detection counter */
      count = 2*nv;
    }
  }

  delete[] V;

  return true;
}

bool Polygon2d::snip(const std::vector<Point2d> &contour, int u, int v, int w, int n, int *V){

  float Ax, Ay, Bx, By, Cx, Cy;

  Ax = contour[V[u]].X();
  Ay = contour[V[u]].Y();

  Bx = contour[V[v]].X();
  By = contour[V[v]].Y();

  Cx = contour[V[w]].X();
  Cy = contour[V[w]].Y();

  // epsilon: 1e-10
  if(1e-10 > (((Bx-Ax)*(Cy-Ay)) - ((By-Ay)*(Cx-Ax))) ) return false;

  for(int p=0;p<n;p++){
    if( (p == u) || (p == v) || (p == w) ) continue;
    if(Triangle2d::inside(contour[V[u]], contour[V[v]], contour[V[w]], contour[V[p]])) return false;
  }

  return true;
}
