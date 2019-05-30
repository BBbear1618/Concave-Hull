# Codes about Alphs shape generation are derived from https://plot.ly/python/alpha-shapes/
# What I did:
#     Find the boundary of an alpha shape incl. inner rings.
# What can be done to improve:
#     Identify rings as outer and inner, and better to return as WKT?

import numpy as np
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.graph_objs import scatter, layout
from plotly import tools as tls
from scipy.spatial import Delaunay

np.seterr(divide='ignore', invalid='ignore')


def sq_norm(v): #squared norm 
    return np.linalg.norm(v)**2

def circumcircle(points,simplex):
    A=[points[simplex[k]] for k in range(3)]
    M=[[1.0]*4]
    M+=[[sq_norm(A[k]), A[k][0], A[k][1], 1.0 ] for k in range(3)]
    M=np.asarray(M, dtype=np.float32)
    S=np.array([0.5*np.linalg.det(M[1:,[0,2,3]]), -0.5*np.linalg.det(M[1:,[0,1,3]])])
    a=np.linalg.det(M[1:, 1:])
    b=np.linalg.det(M[1:, [0,1,2]])
    return S/a,  np.sqrt(b/a+sq_norm(S)/a**2) #center=S/a, radius=np.sqrt(b/a+sq_norm(S)/a**2)
#def distance_2d(p1, p2):
#    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**0.5
#
#def max_length(points, simplex):
#    """
#    return the maximum length of edges of simplex.
#    """
#    return max(distance_2d(points[simplex[0]], points[simplex[1]]),
#               distance_2d(points[simplex[1]], points[simplex[2]]),
#               distance_2d(points[simplex[2]], points[simplex[0]]))


def get_alpha_complex(alpha, points, simplexes):
    #alpha is the parameter for the alpha shape
    #points are given data points 
    #simplexes is the  list of indices in the array of points 
    #that define 2-simplexes in the Delaunay triangulation
    return filter(lambda simplex: circumcircle(points,simplex)[1]<alpha, simplexes)

def Plotly_data(points, complex_s):
    #points are the given data points, 
    #complex_s is the list of indices in the array of points defining 2-simplexes(triangles) 
    #in the simplicial complex to be plotted
    X=[]
    Y=[]
    for s in complex_s:
        X+=[points[s[k]][0] for k in [0,1,2,0]]+[None]
        Y+=[points[s[k]][1] for k in [0,1,2,0]]+[None]
    return X,Y

def Plotly_data_polygon(points, polygon):
    #points are the given data points, 
    #complex_s is the list of indices in the array of points defining 2-simplexes(triangles) 
    #in the simplicial complex to be plotted
    X=[]
    Y=[]
    for loop in polygon:
        for vi in loop:
            X+=[points[vi][0]]
            Y+=[points[vi][1]]
        X+=[None]
        Y+=[None]
    return X,Y

def make_trace(x, y,  point_color='#C0223B', line_color='#404ca0'):# define the trace
                                                                   #for an alpha complex
    return go.Scatter(mode='markers+lines', #set vertices and 
                                         #edges of the alpha-complex
                   name='',
                   x=x,
                   y=y,
                   marker=scatter.Marker(size=6.5, color=point_color),
                   line=scatter.Line(width=1.25, color=line_color),
                  )
    
def make_XAxis(axis_style):
    return layout.XAxis(axis_style)

def make_YAxis(axis_style):
    return layout.YAxis(axis_style)

def plot_delaunay_alphaShape_concaveHull(pts, alpha, file_name):
    """
    Plot Delaunay, Alpha Shape, and Concave Hull in plotly.
    Before using plotly, register first!
    """
    points = pts - np.mean(pts, axis = 0)
    tri = Delaunay(points)
    figure = tls.make_subplots(rows=1, cols=3,
                           subplot_titles=('Delaunay triangulation',  'Alpha shape, alpha={}'.format(alpha), 'Concave hull'),
                           horizontal_spacing=0.1
                          )
    pl_width=1200
    pl_height=460
    title = 'Delaunay triangulation, Alpha Shape, and Concave hull for a Set of 2D Points'
    
    figure['layout'].update(title=title,
                            width=pl_width,
                            height=pl_height
                           )
    
    alpha_complex=get_alpha_complex(alpha, points, tri.simplices)
    
    X1,Y1=Plotly_data(pts, tri.simplices)# get data for Delaunay triangulation
    figure.append_trace(make_trace(X1, Y1), 1, 1)
    
    X2,Y2=Plotly_data(pts, alpha_complex)# data for alpha complex
    figure.append_trace(make_trace(X2, Y2), 1, 2)
    
    X3,Y3=Plotly_data_polygon(pts, concave_hull(pts, alpha))# data for alpha complex
    figure.append_trace(make_trace(X3, Y3), 1, 3)
    
#    for s in alpha_complex: #fill in the triangles of the alpha complex
#        A=pts[s[0]]
#        B=pts[s[1]]
#        C=pts[s[2]]
#        figure['layout']['shapes'].append(dict(path='M '+str(A[0])+',' +str(A[1])+' '+'L '+\
#                                                     str(B[0])+', '+str(B[1])+ ' '+'L '+\
#                                                     str(C[0])+', '+str(C[1])+' Z',
#                                               fillcolor='rgba(173,216,230, 0.5)',
#                                               line=scatter.Line(color='#404ca0', width=1.25),
#                                               xref='x2',
#                                               yref='y2'
#                                               )
#                                         )
    py.plot(figure, filename=file_name, width=1250)
    
def concave_hull(pts, alpha):
    """
    Compute the concave hull of a set of 2D points.
    Input:
        pts: a set of 2D points.
        alpha: threshold for alpha shape generation.
    Output:
        Concave hull of the given points which is represented by a list of several lists of the point indices.
        It may contain inside rings as well as islands.
    """
    points = pts - np.mean(pts, axis = 0)
    tri = Delaunay(points)
#    plt.triplot(points[:, 0], points[:, 1], tri.simplices)
#    print("input pts number: {}, after triangulation: {}".format(len(pts), len(tri.simplices)))
    alpha_complex=list(get_alpha_complex(alpha, points, tri.simplices))
    #record frequency of edges
    f_edge = {}
    for simplex in alpha_complex:
        for j in range(-1, 2):
            vertex1 = simplex[j]
            vertex2 = simplex[j + 1]
            
            # if this edge is found in f_edge,
            # i.e., its duplicate is found in f_edge,
            # then, count + 1
            found = False
            for edge in f_edge.keys():
                v1 = edge[0]
                v2 = edge[1]
                if (vertex1 == v1 and vertex2 == v2) or (vertex1 == v2 and vertex2 == v1):
                    found = True
                    f_edge[edge] += 1
                    break
                    
            # if this edge never appear in f_edge
            if found == False:
                f_edge[(vertex1, vertex2)] = 1
            
#    print(f_edge)
    # remove duplicate edges
    boundary_unordered = [e for e in f_edge.keys() if f_edge[e] == 1]
#    print(boundary_unordered)
    boundary_ordered = []
    oneLoop = []
    v1, v2 = boundary_unordered.pop(0)
    oneLoop.append(v1)
    oneLoop.append(v2)
    while len(boundary_unordered) > 0:
        flag = False
        for e in boundary_unordered:
            # next edge found
            if e[0] == oneLoop[-1]:
                flag = True
                oneLoop.append(e[1])
                boundary_unordered.remove(e)
                break
            
        # next edge not found in left vertices, i.e., one loop formed
        if flag == False:
            boundary_ordered.append(oneLoop)
            oneLoop = []
            v1, v2 = boundary_unordered.pop(0)
            oneLoop.append(v1)
            oneLoop.append(v2)
        elif len(boundary_unordered) == 0:
            boundary_ordered.append(oneLoop)
#    print(boundary_ordered)
    return boundary_ordered
    
    
if __name__ == "__main__":
    filename = "363100012072875_1.txt"
    pts = np.loadtxt(filename)
    plot_delaunay_alphaShape_concaveHull(pts, 0.7, filename.split('.')[0])
#    concave_hull(pts)