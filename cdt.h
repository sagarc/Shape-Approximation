#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <CGAL/basic.h>
#include <CGAL/intersections.h> 
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Object.h>
#include <CGAL/enum.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/centroid.h>


#include "primitives.h"
#include "knn.h"
#include "visilibity/visilibity.hpp"
#include "ramp.h"

#include <stack>
#include<fstream>
#include<utility>



template <class CDT>
class CCDT : public CDT
{
public:

	// typedefs for basic primitives
	typedef typename CDT::Geom_traits Kernel;
	typedef typename Kernel::FT FT;
	typedef typename Kernel::Point_2    Point;
	typedef typename Kernel::Ray_2      Ray;
	typedef typename Kernel::Vector_2   Vector;
	typedef typename Kernel::Segment_2  Segment;
	typedef typename Kernel::Triangle_2 Triangle;
	typedef typename Kernel::Direction_2 Direction;
	typedef typename CGAL::Orientation Orientation;

	typedef typename CCDT<CDT> Cdt;

	// handles
	typedef typename Cdt::Edge    Edge;
	typedef typename Cdt::Face_handle         Face_handle;
	typedef typename Cdt::Vertex_handle       Vertex_handle;

	// iterators
	typedef typename Cdt::Face_iterator       Face_iterator;
	typedef typename Cdt::Edge_iterator       Edge_iterator;
	typedef typename Cdt::Vertex_iterator     Vertex_iterator;

	// circulators
	typedef typename Cdt::Edge_circulator     Edge_circulator;
	typedef typename Cdt::Face_circulator     Face_circulator;
	typedef typename Cdt::Vertex_circulator   Vertex_circulator; 
	typedef typename Cdt::Line_face_circulator   Line_face_circulator; 

	typedef CNeighbor_search<Kernel> KNN;

	typedef typename std::list<Point> Boundary;
	Boundary m_outer_boundary;
	std::list<Boundary> m_inner_boundaries;
	Boundary m_r_outer_boundary;
	std::list<Boundary> m_r_inner_boundaries;
	Boundary m_ptoc_polygon; // only if one hole seen....ptoc-primitive_to_constructed_polyline
	Boundary m_recursive_polygon;
	Boundary inner_rb,outer_b;
	
	Vertex_handle m_max_v,m_curr_v;
	bool m_vis_flag;
	Boundary m_bad_points;
	
	//Boundary m_test_lfc,m_test_lfc_nointersects,m_test_lfc_angle,m_test_lfc_final;
	Point m_next_p;
	//Boundary m_test_vispoly_bad;
	bool m_first_time_flag,m_construction_done;
	double m_total_dist, m_total_dist1;
	
	// all boundary edges
	std::set<Edge> m_all_edges;
	bool m_all_holes_covered;

	//std::list<Point> m_local_maxima;
	//std::list<Point> m_local_maxima_approx;
	//std::list<Face_handle> m_local_maxima_face;
	
	// for dual graph
	typedef std::pair<Vertex_handle,Vertex_handle> V_pair;
	typedef std::pair<Face_handle,Face_handle> F_pair;
	typedef std::pair<F_pair,V_pair> Moat;	
	
	std::list<Segment> m_dual_graph_edges;
	std::list<Point> m_hubs;
	std::list<Face_handle> m_hubfaces;
	std::list<Moat> m_moats;

	
private:
	//Member from Visilibity lib
	typedef typename VisiLibity::Point Vis_point;
	typedef typename VisiLibity::Polygon Vis_polygon;
	typedef typename VisiLibity::Environment Vis_environment;
	typedef typename VisiLibity::Guards Vis_guards;
	typedef typename VisiLibity::Visibility_Polygon V_polygon;
	typedef typename VisiLibity::Visibility_Graph Vis_Graph;
	
	std::vector<Vis_polygon> m_polygons;
	Vis_environment m_environment;
	std::list<Vis_point>m_guards_vec;
	std::list<Vertex_handle> m_internal_vertices;
	Vis_guards m_guards;
	Vis_Graph m_visibility_graph;

public:
	CCDT() 
	{
		m_vis_flag = 0;
	}
	virtual ~CCDT() 
	{
	}

	enum {INSIDE = -1,
		  UNDETERMINED = 0,
		  OUTSIDE = 1,
	      BOUNDARY = 2,
	      HOLE = 3};

public:
	// random (uniform)
	void generators_random(unsigned int nb_generators)
	{
		clear();
		for(unsigned int i=0;i<nb_generators;i++)
			insert(Point(r(),r()));
	}

	// random number between zero and max
	FT r(FT max = 1.0) { return max * (FT)rand() / (FT)RAND_MAX; }

	// OPENGL DRAWINGS

	// draw generators
	void gl_draw_generators(float point_size,
							float red,
							float green,
							float blue)
	{
		::glColor3f(red,green,blue);
		::glPointSize(point_size);
		::glBegin(GL_POINTS);
		Point_iterator p;
		for(p  = points_begin(); 
				p != points_end(); 
				p++) 
			::glVertex2d(p->x(),p->y());
			::glEnd();
	}

	void gl_draw_vertices(float point_size)
	{
		::glPointSize(point_size);
		::glBegin(GL_POINTS);
		Finite_vertices_iterator v;
		for(v  = finite_vertices_begin(); 
				v != finite_vertices_end(); 
				v++) 
		{
			const Point& p = v->point();
			switch(v->location())
			{
			case INSIDE:
				::glColor3ub(0,0,255);
				break;
			case BOUNDARY:
				::glColor3ub(255,0,0);
				break;
			default: // OUTSIDE
				::glColor3ub(0,0,0);
			}
		  
			::glVertex2d(p.x(),p.y());
		}
		::glEnd();
	}

  void gl_draw_unconstrained_edges(const unsigned char r,
                                   const unsigned char g,
                                   const unsigned char b)
  {
	  ::glLineWidth(1.0f);
    ::glColor3ub(r,g,b);

		::glBegin(GL_LINES);
		Finite_edges_iterator e;
		for(e  = finite_edges_begin(); 
				e != finite_edges_end(); 
				e++) 
		{
      Edge edge = *e;
      if(is_constrained(edge))
        continue;
	    Point p1 = edge.first->vertex(ccw(edge.second))->point();
		  Point p2 = edge.first->vertex(cw(edge.second))->point();
			::glVertex2d(p1.x(),p1.y());
			::glVertex2d(p2.x(),p2.y());
		}
		::glEnd();
  }

  void gl_draw_constrained_edges(const unsigned char r,
                                 const unsigned char g,
                                 const unsigned char b)
  {
	::glLineWidth(1.5f);
    ::glColor3ub(r,g,b);
	::glBegin(GL_LINES);
	Finite_edges_iterator e;
	for(e  = finite_edges_begin(); 
			e != finite_edges_end(); 
			e++) 
	{
		Edge edge = *e; // STL pair (see doc)
		if(!is_constrained(edge))
		continue;
		Point p1 = edge.first->vertex(ccw(edge.second))->point();
		Point p2 = edge.first->vertex(cw(edge.second))->point();
		::glVertex2d(p1.x(),p1.y());
		::glVertex2d(p2.x(),p2.y());
	}
	::glEnd();
  }

	// draw inside face
	void gl_draw_face_location(unsigned char r,
		                         unsigned char g,
		                         unsigned char b,
								 const int location)
	{
		::glColor3ub(r,g,b);
		::glBegin(GL_TRIANGLES);
		Finite_faces_iterator f;
		for(f = finite_faces_begin();
			f != finite_faces_end();
			f++)
		{
			if(f->location() == location)
			{
				const Point& p1 = f->vertex(0)->point();
				const Point& p2 = f->vertex(1)->point();
				const Point& p3 = f->vertex(2)->point();
				::glVertex2d(p1.x(),p1.y());
				::glVertex2d(p2.x(),p2.y());
				::glVertex2d(p3.x(),p3.y());
			}
		}
		::glEnd();
	}
	
	void tag_primitives_inside_outside(KNN *pKNN,
		                               const double tolerance)
	{
		tag_faces(pKNN, tolerance);
		tag_vertices_inside_outside();
	}

	void tag_vertices_inside_outside()
	{
		All_vertices_iterator v;
		for(v = all_vertices_begin();
			  v != all_vertices_end();
				v++)
		{
			if(all_locations_around(v,INSIDE))
			v->location() = INSIDE;
			else
			if(all_locations_around(v,OUTSIDE))
				v->location() = OUTSIDE;
			else
			if(all_locations_around(v,HOLE))
				v->location() = OUTSIDE;
			else
				v->location() = BOUNDARY;
		}
	}

	bool all_locations_around(Vertex_handle v, const int location)
	{
		Face_circulator c = incident_faces(v);
		Face_circulator end = c;
		CGAL_For_all(c, end)
			if(c->location() != location)
				return false;
		return true;
	}

	void tag_faces(KNN *pKNN,
		             const double tolerance)
	{
		All_faces_iterator f;
		for(f = all_faces_begin();
			f != all_faces_end();
			f++)
		{
			f->location() = HOLE;
			if(!is_infinite(f))
			{
				const Point& a = f->vertex(0)->point();
				const Point& b = f->vertex(1)->point();
				const Point& c = f->vertex(2)->point();
				Point cc = CGAL::circumcenter(a,b,c);				
				if(inside_tolerance(pKNN, cc, tolerance))
					f->location() = INSIDE;
			}
		}
		tag_face_component_from(infinite_vertex(), HOLE, OUTSIDE);
	}

	void tag_face_component_from(Vertex_handle v,
		                         const int component,
								const int tag)
	{
		Face_handle seed = v->face();
		seed->location() = OUTSIDE;
		std::stack<Face_handle> faces;
		faces.push(seed);

		while(!faces.empty())
		{
			Face_handle f = faces.top();
			faces.pop();
			f->location() = tag;
			const int& location = f->location();
			for(unsigned int i=0;i<3;i++)
				if(f->neighbor(i)->location() == component)
					faces.push(f->neighbor(i));
		}
	}

	bool inside_tolerance(KNN *pKNN,
		                    const Point& query, 
		                    const double tolerance)
	{
		FT d = pKNN->distance_nearest_point(query);
		return d <= tolerance;
	}

	

  // utility function to know if a point is in the "inside" area
  bool inside(const Point& query)
  {
    Face_handle f = locate(query);
    return f->location() == INSIDE;
  }

  // make simple box 
  void add_box(const Point center, const FT size)
  {
    Point a(center.x() - size, center.y() - size);
    Point b(center.x() + size, center.y() - size);
    Point c(center.x() + size, center.y() + size);
    Point d(center.x() - size, center.y() + size);
    insert_constraint(a,b);
    insert_constraint(b,c);
    insert_constraint(c,d);
    insert_constraint(d,a);
  }

	void add_circle(const unsigned int nb_points)
	{
		Point center(0.5,0.5);
		const double incr = 360.0 / double(nb_points);
		std::list<Vertex_handle> vertices;
		for (double angle = 0.0; angle < 360.0; angle += incr)
		{
			double angle_rad = (angle / 360.0) * 6.2831853;
			double x = center.x() + 0.5 * cos(angle_rad);
			double y = center.y() + 0.5 * sin(angle_rad);
			Vertex_handle v = insert(Point(x,y));
			vertices.push_back(v);
		}

		for(int i=0;i<nb_points;i++)
			insert_constraint(vertices[i],vertices[(i+1)%nb_points]);
	}	

	void draw_boundaries()
	{
		// outer
		draw_boundary(m_outer_boundary, 1.0f, 255, 255, 100);

		// inner
		std::list<Boundary>::iterator it;
		for(it = m_inner_boundaries.begin();
			  it != m_inner_boundaries.end();
				it++)
		{
			Boundary& boundary = *it;
			draw_boundary(boundary, 2.0f, 0, 255, 100);		  
		}
	}
	
	void clear_boundaries()
	{
		// outer
		draw_boundary(m_outer_boundary, 2.0f, 0, 0, 0);

		// inner
		std::list<Boundary>::iterator it;
		for(it = m_inner_boundaries.begin();
			  it != m_inner_boundaries.end();
				it++)
		{
			Boundary& boundary = *it;
		  draw_boundary(boundary, 2.0f, 0, 0, 0);
		}
		
	}

	void draw_boundary(Boundary& boundary,
		                 const float line_width,
										 const unsigned char red,
										 const unsigned char green,
										 const unsigned char blue)
	{
	  ::glLineWidth(line_width);
      ::glColor3ub(red, green, blue);

		// draw (closed) edge polyline
		::glBegin(GL_LINE_LOOP);
		std::list<Point>::iterator it;
		for(it = boundary.begin();
			  it != boundary.end();
				it++)
		{
			const Point& p = *it;
			::glVertex2d(p.x(),p.y());
		}
		::glEnd();

		// draw arrows to depict CW vs CCW
		for(it = boundary.begin();
			  it != boundary.end();
				it++)
		{
			const Point& p = *it;

			// pick next point along boundary 
			std::list<Point>::iterator next_it = it;
			next_it++;
			if(next_it != boundary.end())
			{
				const Point& n = *next_it; // next point
				// assemble arrow
				Vector u = n - p;
				Point a = p + 0.7*u;
				Vector u90(-u.y(), u.x());
				Point b = p + 0.3*u + 0.1*u90;
				Point c = p + 0.3*u - 0.1*u90;

				::glBegin(GL_TRIANGLES);
				::glVertex2d(a.x(),a.y());
				::glVertex2d(b.x(),b.y());
				::glVertex2d(c.x(),c.y());
				::glEnd();
			}
		}		
	}

	

	void extract_boundaries(Boundary& m_outer_edge, std::list<Boundary>& m_inner_edges)
	{
		m_all_edges.clear();
		m_all_holes_covered = 0 ;
		extract_outer_boundary(m_outer_edge);
		extract_inner_boundaries(m_inner_edges);
	}

	void extract_inner_boundaries(std::list<Boundary>& m_inner_edges)
	{
		m_inner_boundaries.clear();
		m_inner_edges.clear();
		std::cout << "Extract inner boundaries...";
		Boundary boundary;
		
		Edge seed_edge = find_seed_boundary_edge(HOLE, INSIDE);		
		Edge edge_p(seed_edge);
		while(!m_all_holes_covered)
		{
			boundary.clear();
			do
			{
				boundary.push_back(end_point(edge_p));
				m_all_edges.insert(edge_p);
				edge_p = next_edge_along_interface(edge_p, HOLE, INSIDE);
			}
			while(edge_p != seed_edge);
			m_inner_boundaries.push_back(boundary);
			m_inner_edges.push_back(boundary);
			
			seed_edge = find_seed_boundary_edge(HOLE, INSIDE);	
			edge_p = seed_edge;
		}
		std::cout<<"Extracted "<< m_inner_boundaries.size()<<" inner boundaries\n";
	}
	
	void extract_outer_boundary(Boundary& m_outer_edge)
	{
		std::cout << "Extract outer boundary...";
		m_outer_boundary.clear();
		m_outer_edge.clear();
		Edge seed_edge = find_seed_boundary_edge(OUTSIDE, INSIDE);
		Edge edge = seed_edge;
		do
		{
			m_outer_boundary.push_back(end_point(edge));
			m_outer_edge.push_back(end_point(edge));
			m_all_edges.insert(edge);
			edge = next_edge_along_interface(edge, OUTSIDE, INSIDE);
		}
		while(edge != seed_edge);
		std::cout << "done (" << m_outer_boundary.size() << " points)"<<std::endl;
	}

	Edge next_edge_along_interface(Edge& seed_edge,
		                           const int location, 
								   const int nlocation)
	{
		Edge edge = next(seed_edge);
		while(!is_on_interface_edge(edge,location,nlocation))
		{
			edge = opposite_edge(edge);
			edge = next(edge);
		}
		return edge;
	}

	bool is_on_interface_edge(const Edge& edge, 
		                      const int location, 
							  const int nlocation)
	{
		Face_handle f = edge.first;
		Face_handle nf = (opposite_edge(edge)).first;
		return f->location() == location && nf->location() == nlocation;
	}

	// returns opposite halfedge
	Edge opposite_edge(const Edge& edge)
	{
		Face_handle f = edge.first;
		Face_handle nf = f->neighbor(edge.second);
		Vertex_handle v = f->vertex(ccw(edge.second)); // choose ref. vertex
		int index = nf->index(v);
		return Edge(nf,ccw(index));
	}

	// returns next halfedge with face edge.first
	Edge next(const Edge& edge)
	{
		return Edge(edge.first, cw(edge.second));
	}
	
	Point end_point(const Edge& edge)
    {
		return edge.first->vertex(ccw(edge.second))->point();
	}

	Edge find_seed_boundary_edge(const int location,
		                           const int nlocation)
	{
		All_faces_iterator f;
		for(f = all_faces_begin();
			f != all_faces_end();
			f++)
		{
			if(f->location() != location) // only for faces with specified location
				continue;

			// ckeck for neighboring faces.
			for(int i=0;i<3;i++)
			{
				Face_handle nf = f->neighbor(i);
				if(nf->location() == nlocation && *(m_all_edges.find(Edge(f,i)))==*(m_all_edges.end()))
				{
						return Edge(f,i);
				}			
			}
		}
		m_all_holes_covered  = 1;
		return Edge(Face_handle(),0);
	}



// preprocessing of vertices in environment
//update_environment sets up environment model for visilibity library
//visibility_function() does actual preprocessing of vertices
/*
 * 
 * 
 *All FUNCTIONS FROM HERE ARE WRITTEN BY ME:
 * 
 * 
 * 
 * 
 * */


	void update_environment()
	{
		double epsilon =  0.000000001;
		m_polygons.clear();

	   Vis_polygon polygon;
	   std::list<Point>::iterator it;
		for(it = m_outer_boundary.begin(); it != m_outer_boundary.end(); it++)
		{
			 const Point& p = *it;
			 polygon.push_back(Vis_point(p.x(),p.y()));
		}
		m_polygons.push_back(polygon);

		std::list<Boundary>::iterator bit;
		for(bit = m_inner_boundaries.begin(); bit != m_inner_boundaries.end(); bit++)
		{
			Boundary& boundary = *bit;
			Vis_polygon polygon;
			std::list<Point>::iterator it;
			for(it = boundary.begin(); it != boundary.end(); it++)
			{
			   const Point& p = *it;
			   polygon.push_back(Vis_point(p.x(),p.y()));
			}
			 m_polygons.push_back(polygon);
		}		
		m_environment = Vis_environment(m_polygons);		
	}


	void visibility_function()
	{
		//value set
		double epsilon =  0.000000001;	
		std::cout<< "Visibility function..." << std::endl;	
		outer_b = m_outer_boundary;
		inner_rb = m_inner_boundaries.front();		
		Boundary::iterator outer_nearest_it,inner_nearest_it, vis_it, max_vis_it;		
		
		int count = 0;
		All_vertices_iterator v;
		for(v = all_vertices_begin();
			  v != all_vertices_end();
				v++)
		{				
			if(v->location() == OUTSIDE)
				continue;
			else
			{
				Point& ref_p = v->point();
				Vis_point ref_vp(ref_p.x(),ref_p.y());
				V_polygon v_polygon(ref_vp, m_environment, epsilon);								
				
				std::vector<Vis_point>& vis_poly = v_polygon.vertices_  ;
				std::list<Point> vis_polygon;
				Vis_point vp;Point p;
				for(std::vector<Vis_point>::iterator it = vis_poly.begin(); it != vis_poly.end(); it++)
				{
					vp = *it;
					p = Point(vp.x() , vp.y());
					vis_polygon.push_back(p);
				}
				v->visibility_polygon() = vis_polygon;
				
				outer_nearest_it = find_nearest_point(ref_p,outer_b);
				inner_nearest_it = find_nearest_point(ref_p,inner_rb);
				v->vec_divider() = Vector(*outer_nearest_it , *inner_nearest_it);
				Vector vec_current;
				Point max_vis_p = ref_p;
				for(Boundary ::iterator it = vis_polygon.begin(); it!= vis_polygon.end(); it++)
				{
					vec_current = Vector(ref_p, *it);
					if(CGAL::orientation(v->vec_divider(),vec_current) != CGAL::LEFT_TURN)
					{
						if(CGAL::has_larger_distance_to_point(ref_p, *it , max_vis_p))
						{
							max_vis_p = *it;
						}
					}
				}
				if(max_vis_p == ref_p) m_bad_points.push_back(ref_p);
				v->ccw_max_point() = max_vis_p;
				v->visibility_length() = CGAL::squared_distance(ref_p, max_vis_p);	
				v->visibility_weight() = v->visibility_length();
				v->curr_polyline().push_back(ref_p);
				find_next_candidates(v);
				std::cout<<".";		
			}
		}
		// for greedy approach
		std::cout<<"\n";
		double max_vis_length = 0;
		for(v = all_vertices_begin();
			  v != all_vertices_end();
				v++)
		{
			if(v->location() != OUTSIDE)
			{
				if(v->visibility_length() > max_vis_length)
				{
					max_vis_length = v->visibility_length();
					m_max_v = v;
				}				
			}
		}
		m_vis_flag = 1;
		std::cout<< "done" << std::endl;				
	}	

	void find_next_candidates(Vertex_handle& ref_v)
	{
		//std::cout<<".";						
		Vector vec_curr;	
		Vertex_handle v,v_last;			
		Point& ref_p = ref_v->point();			
		Point& end_p = ref_v->ccw_max_point();	
			
		Point curr_p;
		Line_face_circulator lfc = line_walk(ref_p, end_p);
		Line_face_circulator lfc_end = lfc;
		CGAL_For_all(lfc,lfc_end)
		{			
			if(lfc->location() == INSIDE)
			{
				for(int i=0; i<3; i++)
				{
					v = lfc->vertex(i);
					curr_p = v->point();				

					if(CGAL::angle(end_p,ref_p, curr_p) == CGAL::ACUTE)
					{											
						if( CGAL::bounded_side_2(ref_v->visibility_polygon().begin() ,
							ref_v->visibility_polygon().end(), curr_p) != CGAL::ON_UNBOUNDED_SIDE)								
						{								
							ref_v->next_candidates().push_back(v);
						}
					}					
				}
			}
		}		
	}

	
	void compute_min_nested_polygon()
	{
		std::cout<<"entered min nested\n";
		std::list<Vertex_handle> v_result;	
		int count = 2;
		bool computed = 0;
		while(1)
		{
			computed = 0;
			v_result.clear();
			All_vertices_iterator v;
			for(v = all_vertices_begin();
				  v != all_vertices_end();
					v++)
			{
				if(v->location() != OUTSIDE)
				{
					find_next_vertex_recursive(v);
					if(v->complete() == 1) 
					{
						v_result.push_back(v);
						computed = 1;
					}
				}
			}
			for(v = all_vertices_begin();
				  v != all_vertices_end();
					v++)
			{				
				v->curr_polyline() = v->next_polyline();
			}
			
			if(computed == 1 && v_result.size() > 5)
			{
				std::cout<<"total points are:"<<v_result.size()<<"\n";
				break;
			}
			std::cout<<"\none more loop"<<count<<"\n";
			count++;
		}

		Vertex_handle base_v = v_result.front();		
		m_recursive_polygon = base_v->curr_polyline();
		m_recursive_polygon.push_front(base_v->point());
		
		m_total_dist = 0;
		std::list<Point> :: iterator prev_a,curr_a;
		prev_a = m_recursive_polygon.end();
		prev_a --;
		for( curr_a = m_recursive_polygon.begin(); curr_a!= m_recursive_polygon.end(); curr_a++)
		{
			m_total_dist += CGAL::squared_distance(*curr_a, *prev_a);
			prev_a = curr_a;
		}
		std::cout<<"total length "<<m_total_dist<<"\n";
	}

	
	void find_next_vertex_recursive(Vertex_handle& ref_v)
	{
		std::cout<<".";
		if(ref_v->complete() == 1)
		{
			ref_v->complete() = 2;
			std::cout<<"$";
			return;
		}
		else if(ref_v->complete() == 2);
		else
		{		
			bool max_assigned = 0;
			double max_weight = 0,curr_weight;		
			Vector vec_curr;	
			Vertex_handle max_vertex,v,v_last;	
			max_vertex = ref_v;
			Point& ref_p = ref_v->point();			
			Point& end_p = ref_v->ccw_max_point();		
			
			Point curr_p, curr_last_p;
			for(std::list<Vertex_handle> :: iterator it = ref_v->next_candidates().begin();
				it != ref_v->next_candidates().end(); it++)
			{
				v = *it;
				curr_p = v->point();
				curr_last_p = (v->curr_polyline()).back();						
				curr_weight = CGAL::squared_distance(ref_p, curr_p) + v->visibility_weight();
				if(curr_weight > max_weight)	
				{							
					vec_curr = Vector(ref_p, curr_last_p);
					if( ref_v->intersect_flag() &&
						CGAL::orientation(ref_v->vec_divider(), vec_curr) != CGAL::RIGHT_TURN &&
						CGAL::bounded_side_2(ref_v->visibility_polygon().begin(), 
											ref_v->visibility_polygon().end(), curr_last_p) != CGAL::ON_UNBOUNDED_SIDE)
					{									
						ref_v->next_polyline() = v->curr_polyline();
						(ref_v->next_polyline()).push_front(curr_p);
						ref_v->complete() = 1;								
						max_assigned = 1;
						return;
					}															
					max_weight = curr_weight;									
					max_vertex =  v;
					max_assigned = 1;
				}			
			}
			//std::cout<<"over\n";
			if(max_vertex != ref_v)
			{
				ref_v->visibility_weight() = max_weight;				
				ref_v->next_polyline() = max_vertex->curr_polyline();
				ref_v->next_polyline().push_front(max_vertex->point());

				if(!ref_v->intersect_flag())
				{
					if(intersects_boundary(Segment(ref_p, (ref_v->next_polyline()).back()) , inner_rb))
					{	
						ref_v->intersect_flag() = 1;						
					}
				}
			}
			else //!max_assigned
			{				
				std::cout<<"error ";				
			}
		}
	}	

	//Finds next potential point as vertex of polygon
	Vertex_handle find_next_vertex(Vertex_handle& ref_v)
	{
		bool max_assigned = 0;
		double max_weight = 0,curr_weight;		
		Vector vec_current;
		Vertex_handle max_vertex,v;
		
		Point max_p = m_max_v->point();
		Point ref_p = ref_v->point();
		Point curr_p;
		Point& end_p = ref_v->ccw_max_point();
		Line_face_circulator lfc = line_walk(ref_p, end_p);
		Line_face_circulator lfc_end = lfc;
		CGAL_For_all(lfc,lfc_end)
		{			
			if(lfc->location() == INSIDE)
			{
				for(int i = 0;i<3;i++)
				{
					v = lfc->vertex(i);
					curr_p = v->point();						
					if(CGAL::angle(end_p,ref_p, curr_p) == CGAL::ACUTE)
					{							
						if( CGAL::bounded_side_2(ref_v->visibility_polygon().begin() ,
													ref_v->visibility_polygon().end(), curr_p) != CGAL::ON_UNBOUNDED_SIDE)								
						{								
							vec_current = Vector(curr_p, m_max_v->point());
							if( !m_first_time_flag &&
								CGAL::orientation(v->vec_divider(), vec_current) != CGAL::LEFT_TURN &&
								CGAL::bounded_side_2(v->visibility_polygon().begin() ,
													v->visibility_polygon().end(), max_p) != CGAL::ON_UNBOUNDED_SIDE)
							{									
								max_vertex = v;
								max_assigned = 1;
								return max_vertex;
							}
							curr_weight = CGAL::squared_distance(ref_p, curr_p) + v->visibility_length();
							if(curr_weight > max_weight)
							{
								max_weight = curr_weight;									
								max_vertex =  v;
								max_assigned = 1;
							}
						}
					}				
				}
			}
		}	

		if(max_assigned)
			return max_vertex;
		else
		{	
			std::cout<<"error - no max_vertex assigned\n";
			std::cout<<ref_p << " "<<end_p<<"\n";							
		}
	}
		
	void construct_polygon()
	{		
		int count = 1;
		start_construct_polygon();
		while(!m_construction_done)
		{
			construct_polygon_next_point();
			count++;
		}
		m_total_dist1 = 0;
		std::list<Point> :: iterator prev_a,curr_a;
		prev_a = m_ptoc_polygon.end();
		prev_a--;
		for( curr_a = m_ptoc_polygon.begin(); curr_a!= m_ptoc_polygon.end(); curr_a++)
		{
			m_total_dist1 += CGAL::squared_distance(*curr_a, *prev_a);
			prev_a = curr_a;
		}
		std::cout<<"total length in ptoc is"<<m_total_dist1<<"\n";
		std::cout<<"no of segments are: "<< count<<"\n";

	} 

	void start_construct_polygon()
	{
		//first point  is known to be m_max_v
		std::cout<<"started coonstruct_polygon\n";
		m_first_time_flag = 1;
		m_construction_done = 0;		
		Vertex_handle curr_v,next_v;
		Point next_p,curr_p,max_p;
		m_curr_v = m_max_v;		
		next_v = find_next_vertex(m_curr_v);
		
		max_p = m_max_v->point();	
		curr_p = m_curr_v->point();
		next_p = next_v->point();

		std::cout<<"added point "<<curr_p<<"\n";
		m_ptoc_polygon.push_back(curr_p);
		m_curr_v = next_v;		
		curr_p = m_curr_v->point();
		std::cout<<"added point "<<curr_p<<"\n";
		m_ptoc_polygon.push_back(m_curr_v->point());
	}			

	void construct_polygon_next_point()
	{
		std::cout<<"entered next coonstruct_polygon\n";		
		Vertex_handle next_v;
		Point curr_p, max_p;
		max_p = m_max_v->point();
		curr_p = m_curr_v->point();

		if(intersects_boundary(Segment(curr_p,max_p),inner_rb))
			m_first_time_flag = 0;		
		
		Vector vec_current = Vector(curr_p, m_max_v->point());
		Vector vec_divider = m_max_v->vec_divider();			
		if(CGAL::orientation(vec_divider,vec_current) != CGAL::LEFT_TURN &&
			CGAL::bounded_side_2(m_curr_v->visibility_polygon().begin(),
								m_curr_v->visibility_polygon().end(), max_p) != CGAL::ON_UNBOUNDED_SIDE)			
		{				
			m_construction_done = 1;
			std::cout<<"done constructing\n";
		}
		else
		{			
			next_v = find_next_vertex(m_curr_v);
			m_ptoc_polygon.push_back(next_v->point());
			m_curr_v = next_v;
			std::cout<<"added point "<<m_curr_v->point()<<"\n";
		}			
	}	


	Point find_circumcenter(Face_handle& f)
	{
		Point& a = f->vertex(0)->point();
		Point& b = f->vertex(1)->point();
		Point& c = f->vertex(2)->point();
		Point cc = CGAL::circumcenter(a,b,c);
		return cc;
	}
	
	void find_dual_graph()
	{
		All_faces_iterator f;
		for(f = all_faces_begin();
			f != all_faces_end();
			f++)
		{
			if(f->location() == INSIDE)
			{
				Point cc = find_circumcenter(f);				
				int deg_count = 0;
				for(int i=0;i<3;i++)
				{
					Face_handle nf = f->neighbor(i);
					if(nf->location() == INSIDE)
					{
						Point cc1 =find_circumcenter(nf);
						m_dual_graph_edges.push_back(Segment(cc,cc1));
						deg_count++;					
					}
				}
				f->degree() = deg_count;
			}
		}

		// to remove degre one faces
		bool still_flag= 0;
		do
		{
			still_flag = 0;
			for(f = all_faces_begin();
				f != all_faces_end();
				f++)
			{
				if(f->location() == INSIDE && f->degree() == 1)
				{				
					f->degree()--;
					still_flag = 1;
					for(int i=0;i<3;i++)
					{
						Face_handle nf = f->neighbor(i);
						if(nf->location() == INSIDE && nf->degree()!=0)
						{	
							nf->degree()--;
						}
					}
				}
			}
		}
		while(still_flag == 1);

		// finding of moats and cycle of degree 2 vertices!!
		for(f = all_faces_begin();
				f != all_faces_end();
				f++)
		{
			if(f->location() == INSIDE && f->degree() == 3)
			{
				Point cc = find_circumcenter(f);
				m_hubs.push_back(cc);
				f->moat_count() = 3;
				m_hubfaces.push_back(f);
			}			
		}
		std::cout<<"\n";

		

		std::list<Face_handle> ::iterator f_handle;
		int moat_visit_no = 0;
		for(f_handle = m_hubfaces.begin(); f_handle!= m_hubfaces.end(); f_handle++)
		{
			Face_handle face = *f_handle;			
			if(face->moat_count()!=0)
			{
				for(int i=0;i<3;i++)
				{
					moat_visit_no++;
					if(face->neigh(i) == 0)
					{
						Vertex_handle va = face->vertex(cw(i));
						Vertex_handle vb = face->vertex(ccw(i));
					
						Face_handle t_face;
						Face_handle p_face = face;
						Face_handle n_face = face->neighbor(i);
					    n_face->visit_no() = moat_visit_no;
						while(n_face->degree() < 3)
						{
							bool check = 0;
							for(int j=0;j<3;j++)
							{
								t_face = n_face->neighbor(j);
								if(t_face->location() == INSIDE && 
									t_face->degree()>=2 && 
									t_face->visit_no()!= moat_visit_no &&
									t_face != n_face)
								{
									p_face = n_face;
									n_face = t_face;
									n_face->visit_no() = moat_visit_no;
									std::cout<<find_circumcenter(n_face)<<" ";
									check = 1;
									break;
								}							
							}
							if(!check) std::cout<<"$";
						}
						
						std::cout<<"\n";
						
						for(int j=0;j<3;j++)
						{
							if(n_face->neighbor(j) == p_face || n_face->visit_no() == moat_visit_no )
							{
								std::cout<<"marked end hub\n";
								n_face->neigh(j) = 1;
								n_face->moat_count()--;
								break;
							}
						}
						face->neigh(i) = 1;
						face->moat_count()--;

						//first face in F_pair is source_face and second is sink face
						//first vertex in V_pair corresponds to outer_boundary and second corresponds to inner_boundary
						m_moats.push_back(Moat(F_pair(face,n_face),V_pair(va,vb)));					
					}
				}
			}
		}
		std::cout<<"total hubs are:"<<m_hubs.size()<<"\n";
		std::cout<<"total hubsfaces are:"<<m_hubfaces.size()<<"\n";
		std::cout<<"total moats are"<<m_moats.size()<<"\n";
		std::cout<<"total segment in graph are:"<<m_dual_graph_edges.size()<<"\n";
		
		std::list<Moat> ::iterator moat_it;
		for(moat_it = m_moats.begin(); moat_it!= m_moats.end();moat_it++)
		{
			Moat t_moat = *moat_it;
			V_pair t_v_pair = t_moat.second;
			F_pair t_f_pair = t_moat.first;
			Face_handle f_handle = t_f_pair.first;			
			Point cc = find_circumcenter(f_handle);
			Point ob = (t_v_pair.first)->point();
			Point ib = (t_v_pair.second)->point();
			f_handle = t_f_pair.second;			
			Point cc1 = find_circumcenter(f_handle);
			std::cout<<"moat from " << cc <<" to " << cc1 <<"\n"<<"bank along"<<ob<<" "<<ib<<"\n";
			
		}
	}


	void process_moat(Moat& moat)
	{
		F_pair f_pair = moat.first;
		V_pair v_pair = moat.second;
		Point outer_nearest_v = v_pair.first;
		Point inner_nearest_v = v_pair.second;
		Face_handle source_moat = f_pair.first;
		face_handle target_moat = f_pair.second;
		Point& a = f->vertex(0)->point();
		Point& b = f->vertex(1)->point();
		Point& c = f->vertex(2)->point();
		std::list<Point> source_tr;
		source_tr.push_back(a);
		source_tr.push_back(b);
		source_tr.push_back(c);

		double max_area;
		Vertex_handle source_v;
		All_vertices_iterator v;
		for(v = all_vertices_begin();
			  v != all_vertices_end();
				v++)
		{
			if( CGAL::bounded_side_2(source_tr.begin() ,
							source_tr.end(), v->point()) != CGAL::ON_UNBOUNDED_SIDE)
			{
				Point curr_p  = v->point();
				Vis_point curr_vp(curr_p.x(),curr_p.y());
				V_polygon v_polygon(curr_vp, m_environment, epsilon);
				std::vector<Vis_point>& vis_poly = v_polygon.vertices_  ;
				std::list<Point> vis_polygon;
				Vis_point vp;Point p;
				for(std::vector<Vis_point>::iterator it = vis_poly.begin(); it != vis_poly.end(); it++)
				{
					vp = *it;
					p = Point(vp.x() , vp.y());
					vis_polygon.push_back(p);
				}
				v->visibility_polygon() = vis_polygon;
				v->visibility_area() = v_polygon.area();
				if(v->visibility_area() > max_area) 
				{
					max_area = v->visibility_area() ;
					source_v = v;
				}
			}
		}

		Point ref_p = source_v->point();		
		source_v->vec_divider() = Vector(outer_nearest_v->point(), inner_nearest_v->point());
		Vector vec_current;
		Point max_vis_p = ref_p;
		for(Boundary ::iterator it = (source_v->visibility_polygon()).begin(); 
			it!= (source_v->visibility_polygon()).end(); it++)
		{
			vec_current = Vector(ref_p, *it);
			if(CGAL::orientation(source_v->vec_divider(),vec_current) != CGAL::LEFT_TURN)
			{
				if(CGAL::has_larger_distance_to_point(ref_p, *it , max_vis_p))
				{
					max_vis_p = *it;
				}
			}
		}
		if(max_vis_p == ref_p)std::cout<<"error in finding max_vis_length";
		source_v->ccw_max_point() = max_vis_p;
		v->visibility_length() = CGAL::squared_distance(ref_p, max_vis_p);	
		v->visibility_weight() = v->visibility_length();
		v->curr_polyline().push_back(ref_p);
	}



	
	double my_area(Vis_point& vp)
	{
		double epsilon = 0.00000001;
		V_polygon v_polygon(vp, m_environment, epsilon);
		return v_polygon.area();
	}

	void gl_draw_vis_area(Ramp& ramp)
	{
		// compute max viz area to calibrate the ramp
		double max_length = 0.0;
		All_vertices_iterator v;
		for(v = all_vertices_begin();
			  v != all_vertices_end();
				v++)
		{
			if(v->location() == OUTSIDE)
				continue;
			const double length = v->visibility_length();
			max_length = length > max_length ? length : max_length;
		}
		if(max_length == 0.0)
			return;

    ::glBegin(GL_TRIANGLES);
		Finite_faces_iterator f;
		for(f = finite_faces_begin();
			f != finite_faces_end();
			++f)
		{
			if(f->location() == INSIDE)
			{
				const Vertex_handle& v0 = f->vertex(0);
				const Vertex_handle& v1 = f->vertex(1);
				const Vertex_handle& v2 = f->vertex(2);

				const Point& p0 = v0->point();
				const Point& p1 = v1->point();
				const Point& p2 = v2->point();
				
				ramp.gl_color(v0->visibility_length(), max_length);
				::glVertex2d(p0.x(),p0.y());

				ramp.gl_color(v1->visibility_length(), max_length);
				::glVertex2d(p1.x(),p1.y());

				ramp.gl_color(v2->visibility_length(), max_length);
				::glVertex2d(p2.x(),p2.y());
			}
		}
		::glEnd();
	}

	

	bool intersects_boundary(Segment& s1 , Boundary& b)
	{
		Boundary :: iterator  it;
		Point prev,curr,ipoint;
		if(b.size() <=1) return false; 
		it = b.begin(); curr = *it;
		it++;
		Segment s_curr;
		for(it; it!=b.end(); it++)
		{
			prev = curr;
			curr = *it;				
			if(prev == s1.source() || curr == s1.source() || prev == s1.target() || curr == s1.target()) continue;
			else
			{
				s_curr = Segment(prev,curr);
				if(do_intersect(s_curr,s1)) return true;
			}			
		}
		return false;
	}

	typename Boundary::iterator find_nearest_point(Point& base_point , Boundary& b)
	{
		Boundary :: iterator it,it_min;
		it = b.begin();it_min = it;
		Point curr_min = *it,curr; 		
		for(it = b.begin(); it!=b.end(); it++)
		{
			curr =*it;
			if(CGAL::has_smaller_distance_to_point(base_point, curr, curr_min))
			{	
				curr_min = curr;
				it_min = it;
			}
		}
		return it_min;
	}
	
	
	void draw_primitive()
	{
		::glLineWidth(0.2f);
		::glColor3ub(0, 0, 255);			
		::glBegin(GL_LINE_LOOP);
		std::list<Point>:: iterator it;
		for(it = m_ptoc_polygon.begin(); it != m_ptoc_polygon.end(); it++)
		{			
			const Point& p = *it;
			::glVertex2d(p.x(),p.y());
		}
		::glEnd();
		::glColor3ub(0, 255, 0);
		::glBegin(GL_LINE_LOOP);
		for(it = m_recursive_polygon.begin(); it != m_recursive_polygon.end(); it++)
		{			
			const Point& p = *it;
			::glVertex2d(p.x(),p.y());
		}
		::glEnd();	

		if(m_vis_flag)
		{
			::glColor3ub(120,120,220);
			::glPointSize(5.0f);
			::glBegin(GL_POINTS);			
			::glColor3ub(220,120,220);
			for(Boundary::iterator it = m_bad_points.begin(); it!= m_bad_points.end(); it++)
			{
				::glVertex2d((*it).x(),(*it).y());
			}

			::glEnd();			
		}
		

	
		
		::glLineWidth(0.2f);
		::glColor3ub(120, 0, 120);	
		Point p1,p2;
		for(std::list<Segment> :: iterator it = m_dual_graph_edges.begin(); it != m_dual_graph_edges.end(); it++)
		{
			Segment& s= *it;					
			::glBegin(GL_LINE_LOOP);
			p1 = s.source();
			::glVertex2d(p1.x(),p1.y());
			p2 = s.target();
			::glVertex2d(p2.x(),p2.y());
			::glEnd();	
		}

		::glColor3ub(220, 120, 220);	
		::glPointSize(6.0f);
		::glBegin(GL_POINTS);		
		for(Boundary::iterator it = m_hubs.begin(); it!= m_hubs.end(); it++)
		{
			::glVertex2d((*it).x(),(*it).y());
		}		
		::glEnd();	

		std::list<Moat> ::iterator moat_it;
		for(moat_it = m_moats.begin(); moat_it!= m_moats.end();moat_it++)
		{
			Moat t_moat = *moat_it;
			V_pair t_v_pair = t_moat.second;
			F_pair t_f_pair = t_moat.first;
			Face_handle f_handle = t_f_pair.first;			
			Point cc = find_circumcenter(f_handle);
			Point ob = (t_v_pair.first)->point();
			Point ib = (t_v_pair.second)->point();
			::glLineWidth(0.2f);
			::glColor3ub(0,120, 250);
			::glBegin(GL_LINE_LOOP);
			::glVertex2d(cc.x(),cc.y());
			::glVertex2d(ob.x(),ob.y());
			::glVertex2d(ib.x(),ib.y());
			::glEnd();		
			
		}		
	}
};

#endif // TRIANGULATION_H
