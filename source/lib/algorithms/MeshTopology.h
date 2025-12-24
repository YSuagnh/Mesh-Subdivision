#ifndef MESH_TOPOLOGY_H
#define MESH_TOPOLOGY_H

#include "tiny_obj_loader.h"
#include <deque>
#include <array>
#include <iterator>
#include <math.h>

#if defined(__GNUG__)
    #define ATTRIBUTE_NOINLINE [[gnu::noinline]]
#elif defined(__clang__)
    #define ATTRIBUTE_NOINLINE [[clang::noinline]]
#elif defined(_MSC_VER)
    #define ATTRIBUTE_NOINLINE [[msvc::noinline]]
#else
    #define ATTRIBUTE_NOINLINE 
#endif

/*
Possible references:
1. **mesh topology** and Loop subdivision: https://pbr-book.org/3ed-2018/Shapes/Subdivision_Surfaces
2. root3 subdivision: https://dl.acm.org/doi/10.1145/344779.344835
3. mesh optimization: https://sites.stat.washington.edu/wxs/Siggraph-93/siggraph93.pdf
4. other. "Polygon Mesh Processing": http://staff.ustc.edu.cn/~lgliu/Courses/DGP_2014_autumn-winter/References/Book_Polygon%20Mesh%20Processing.pdf

Some known constrains apply to the mesh:
1. Every edge in mesh should have 2 or less neighbor faces.
    If an edge in mesh got 3 or more neighbor faces, the program will report "Mesh not valid"
    and exit.
2. Mesh should be consistently ordered. Mobius strip cannot be consistently ordered,
    so it would not be supported (check chapter "Mesh Representation" in webpage:
    https://pbr-book.org/3ed-2018/Shapes/Subdivision_Surfaces).
    In the program I presume this situation will never happen.
3. There should be no holes in surface. For example:
                                    /\
                                   /  \
                                  /    \
                                 /face1 \
                             v_2/________\ v_1
                               /\        /\      
                              /  \ hole /  \   
                             /    \    /    \
                            /face2 \  /face3 \ 
                           /________\/________\
                                    v_3     
    For such a mesh, we cannot find v_1's all neighbor faces(face1, face3) by FindFsByVtx().
    Similarly, v_1's neighbor vertices also cannot be found with FindVtxNeighbors().
    Some models downloaded on internet will have such a structure.
*/              

#define NEXT(i) ((i + 1) % 3)
#define PREV(i) ((i + 2) % 3)


namespace meshTopology
{
typedef tinyobj::real_t real_t;
typedef std::vector<real_t>::size_type r_size_type;

struct SDFace;
struct SDEdge;
struct Point;
struct SDAttrib;
struct SDShape;
struct SDVertex;
struct SDNormal;
struct SDTexcoord;

struct Point
{
    std::array<real_t, 3> p{0};

    Point(real_t a, real_t b, real_t c): p{a, b, c}{}
    Point(): Point(0, 0, 0){}

    real_t &operator[] (const int &i) {return p[i];}
    const real_t &operator[] (const int &i) const {return p[i];}
    real_t operator* (const Point &a) const
    {
        return a[0]*p[0] + a[1]*p[1] + a[2]*p[2];
    }
    Point operator* (const real_t &a) const
    {
        return Point(a*p[0], a*p[1], a*p[2]);
    }
    Point operator/ (const real_t &a) const
    {
        return Point(p[0]/a, p[1]/a, p[2]/a);
    }
    Point operator+ (const Point &a) const
    {
        return Point(a[0]+p[0], a[1]+p[1], a[2]+p[2]);
    }
    Point operator- (const Point &a) const
    {
        return Point(p[0]- a[0], p[1] - a[1], p[2] - a[2]);
    }
    bool operator== (const Point &a) const
    {
        return p[0] == a.p[0] && p[1] == a.p[1] && p[2] == a.p[2];
    }
};
Point operator* (const real_t &t, const Point &p)
{
    return Point(t*p[0], t*p[1], t*p[2]);
}
Point absPoint(const Point &a)
{
    return Point(fabs(a[0]), fabs(a[1]), fabs(a[2]));
}
std::ostream &operator<< (std::ostream &os, const Point &p)
{
    std::ostream_iterator<real_t> outIter(os, " ");
    os << "Point{";
    std::copy(p.p.begin(), p.p.end(), outIter);
    os << "}";
    return os;
}
ATTRIBUTE_NOINLINE void pInfo(const Point&);
void pInfo(const Point &p)
{
    std::cout << p << std::endl;
}

struct SDVertex 
{
    int valence = 0;
    int index = 0;
    SDVertex *child = nullptr;
    SDFace *startFace = nullptr;
    Point p;             // vertex indices #TODO: what about texcoord index/normal index?
    SDVertex(real_t a, real_t b, real_t c): p{a, b, c}{}
    SDVertex(const Point &a) :p{a} {}
    SDVertex() = default;

    real_t &operator[] (int i) {return p[i];}
    real_t const& operator[] (int i) const {return p[i];}
};

struct SDNormal
{
    Point p;
    int index = 0;
    SDNormal(real_t a, real_t b, real_t c): p{a,b,c}{}

    real_t &operator[] (const int &i) {return p[i];}
    real_t const& operator[] (const int &i) const {return p[i];}
    real_t operator* (const SDNormal &n)
    {
        return (p * n.p);
    }
};

struct SDTexcoord
{
    Point p;
    int index = 0;
    SDTexcoord(real_t a, real_t b, real_t c): p{a,b,c}{}
};

struct SDAttrib
{
    std::vector<SDVertex*> vs{};
    std::vector<SDNormal*> ns{};
    std::vector<SDTexcoord*> ts{};
};
std::ostream &operator<< (std::ostream &os, const SDAttrib &a)
{
    std::ostream_iterator<Point> outIter(os, "\n");
    os << "Attrib {\n";
    std::for_each(a.vs.begin(), a.vs.end(), [&](SDVertex *v){*outIter++=v->p;});
    os << "};";
    return os;
}

ATTRIBUTE_NOINLINE void pInfo(const SDAttrib&);
void pInfo(const SDAttrib &a)
{
    std::cout << a << std::endl;
}


struct SDEdge 
{
    /*
    * NOTE: 1.For boundary edge, fs[1] is NULL. 
    *   For any legal edge, fs[0] will not be NULL 
    * 2.vertex in verts may NOT pretain the order in which 
    *   they are passed in constructor
    */
    // primary attributes
    std::array<SDFace*, 2> fs{0};   
    std::array<SDVertex*, 2> verts{0};
    // attributes for mesh algorithms
    int f0edgeNum = 0, f1edgeNum = 0;
    int index{0};
    SDEdge(SDVertex *v0, SDVertex *v1)
    {
        verts[0]=std::min(v0, v1); 
        verts[1]=std::max(v0, v1);
    }
    inline bool boundary() const {return fs[1]? false: true;}
    SDVertex *operator[] (const int i) {return verts[i];}
    const SDVertex *operator[] (const int i) const{return verts[i];}
};
bool operator< (const SDEdge &e1, const SDEdge &e2)
{
    if (e1.verts[0] == e2.verts[0]) return (e1.verts[1] < e2.verts[1]);
    return (e1.verts[0] < e2.verts[0]);
}
bool operator== (const SDEdge &e1, const SDEdge &e2)
{
    return (e1.verts[0]==e2.verts[0]) && (e1.verts[1]==e2.verts[1]);
}
std::ostream &operator<< (std::ostream &os, SDEdge e)
{
    os << "Edge {";
    os << e[0]->p << ", " << e[1]->p;
    os << "};";
    return os;
}
ATTRIBUTE_NOINLINE void pInfo(const SDEdge&);
void pInfo(const SDEdge &e)
{
    std::cout << e << std::endl;
}

struct SDShape
{
    bool smooth = false;
    std::vector<SDFace*> fs{};
    std::set<SDEdge> es{};
    ~SDShape()
    {
    }
};
std::ostream &operator<< (std::ostream &os, const SDShape &s)
{
    std::ostream_iterator<SDFace> outIter(os, "\n");
    std::ostream_iterator<SDEdge> outIter1(os, "\n");
    os << "Shape {\n";
    std::for_each(s.fs.begin(), s.fs.end(), [&](SDFace *f){*outIter++=*f;});
    std::for_each(s.es.begin(), s.es.end(), [&](SDEdge e){*outIter1++=e;});
    os << "};";
    return os;
}
ATTRIBUTE_NOINLINE void pInfo(const SDShape&);
void pInfo(const SDShape &s)
{
    std::cout << s << std::endl;
}

struct SDFace 
{
    // primary attributes
    std::array<SDVertex*, 3> verts{0};     // even vertices belong to face
    std::array<SDNormal*, 3> vns{0};
    std::array<SDTexcoord*, 3> vts{0};
    std::array<SDFace*, 3> neighborFs{0};          // neighbor faces

    // attributes for mesh algorithms
    std::array<SDVertex*, 3> sdverts{0};   // odd vertices created by subdivision
    std::array<SDNormal*, 3> sdvns{0};     // normals for odd vertices
    std::array<SDTexcoord*, 3> sdvts{0};   // texcoords for odd vertces

    std::array<SDFace*, 4> child{0};    // subdivision faces
    int index{-1};
    SDFace(SDVertex *a, SDVertex *b, SDVertex *c): verts{a, b, c}{}
    SDFace() = delete;
    SDFace *prevFace(const SDVertex *v) {return neighborFs[PREV(getVertIndex(v))];}
    SDFace *nextFace(const SDVertex *v) {return neighborFs[getVertIndex(v)];}
    inline int getVertIndex(const SDVertex *v)
    {
        for (int i = 0; i < 3; i++)
            if (verts[i] == v)
                return i;
        assert(0);  // raise an error if v is not found in this->verts
        return -1;
    }
    inline SDVertex *otherVertex(const SDEdge *e)
    {
        for (int i = 0; i < 3; i++)
            if (verts[i] != e->verts[0] && verts[i] != e->verts[1])
                return verts[i];
        return nullptr;     // this should never be reached
    }
    inline SDFace *otherNeighbor(const SDVertex *v)
    {
        return neighborFs[NEXT(getVertIndex(v))];
    }
    SDVertex *operator[] (int i) {return verts[i];}
    const SDVertex *operator[] (int i) const{return verts[i];}

    ~SDFace()
    {
    }
};
std::ostream &operator<< (std::ostream &os, const SDFace fs)
{
    os << "Face {P: ";
    os << fs[0]->p << "," << fs[1]->p << "," << fs[2]->p;
    os << "};";
    return os;
}

ATTRIBUTE_NOINLINE void pInfo(const SDFace&);
void pInfo(const SDFace &f)
{
    std::cout << f << std::endl;
}

bool isBoundary(SDVertex *v)
{
    SDFace *f = v->startFace;
    // std::cout << "startface: " << *f << std::endl;
    do {
        // std::cout << *f << std::endl;
        f = f->nextFace(v);
    } while (f && f != v->startFace);
    if (f)
        return false;
    else
        return true;
}

bool isBoundary(SDEdge e)
{
    return e.fs[1] ?true: false;
}

/*
* input: v
* output: range begining at d_first
* describe:
*   find all faces surrounding the vertex v. after process, 
*   faces will be output to range begining at d_first. 
*   d_first should meet requirements of LegacyOutputIterator
*/
template<class OutputIt>
static void findFsByVtx(const SDVertex *v, OutputIt d_first)
{
    std::deque<SDFace*> fs;
    SDFace *cur = v->startFace;
    do {
        fs.push_back(cur);
        cur = cur->nextFace(v);
    } while (cur && cur != v->startFace);

    if (!cur)
    {
        cur = v->startFace;
        cur = cur->prevFace(v);
        while (cur && cur != v->startFace)
        {
            fs.push_front(cur);
            cur = cur->prevFace(v);
        }
    }
    std::copy(fs.begin(), fs.end(), d_first);
}

/*
Find neighbor vertices of a given vertex v.
Similar to findFsByVtx().
*/
template<class OutputIt>
static void findVtxNeighbors(const SDVertex *v, OutputIt d_first)
{
    std::deque<SDVertex*> neighbors;
    SDFace *cur = v->startFace;
    do {
        neighbors.push_back((*cur)[NEXT(cur->getVertIndex(v))]);
        cur = cur->nextFace(v);
    } while (cur && cur != v->startFace);

    if (!cur)
    {
        cur = v->startFace;
        do {
            neighbors.push_front((*cur)[PREV(cur->getVertIndex(v))]);
            cur = cur->prevFace(v);
        } while (cur);
    }
    std::copy(neighbors.begin(), neighbors.end(), d_first);
}

/*
* attrib: input
* shapes: input
* sdattrib: output
* sdshapes: output
*/
static void createSDShapes(const tinyobj::attrib_t &attrib, const std::vector<tinyobj::shape_t> &shapes,
    SDAttrib &sdattrib, std::vector<SDShape> &sdshapes)
{
    std::vector<SDVertex*> verts;
    std::vector<SDNormal*> vns;
    std::vector<SDTexcoord*> vts;

    // build data structures from input
    if (attrib.vertices.size() % 3 != 0)
        std::cout << "ERROR: loopSubdivision(): wrong size\n";
    for (std::vector<real_t>::size_type i = 0; i < attrib.vertices.size() / 3; i++)
    {
        verts.emplace_back(new SDVertex(attrib.vertices[3*i], attrib.vertices[3*i+1], attrib.vertices[3*i+2]));
        verts.back()->index = verts.size() - 1;
    }
    for (std::vector<real_t>::size_type i = 0; i < attrib.normals.size() / 3; i++)
        vns.emplace_back(new SDNormal(attrib.normals[3*i], attrib.normals[3*i+1], attrib.normals[3*i+2]));
    for (std::vector<real_t>::size_type i = 0; i < attrib.texcoords.size() / 3; i++)
        vts.emplace_back(new SDTexcoord(attrib.texcoords[3*i], attrib.texcoords[3*i+1], attrib.texcoords[3*i+2]));

    for (std::vector<tinyobj::shape_t>::size_type shapei = 0; shapei < shapes.size(); shapei++)
    {
        sdshapes.emplace_back();
        std::vector<SDFace*> &fs = sdshapes[shapei].fs;
        std::set<SDEdge> edges;
        // convert to SD data structure
        const tinyobj::shape_t &shape = shapes[shapei];
        const std::vector<unsigned int> &numFaceVerts = shape.mesh.num_face_vertices;
        const std::vector<tinyobj::index_t> &indices = shape.mesh.indices;
#ifndef NDEBUG
        std::cout << "numfaceverts.size: " << numFaceVerts.size() << "\n";
        std::cout << "verts.size: " << verts.size() << "\n";
        std::cout << "indices.size: " << indices.size() << "\n";
#endif
        for (std::vector<unsigned int>::size_type i = 0; i < numFaceVerts.size(); i++)
        {
            static int a = 0;
            if (numFaceVerts[i] != 3)
            {
                std::cout << "ERROR: loopSubdivision(): need triangle mesh/n";
                return;
            }

            fs.push_back(new SDFace(verts[indices[3*i].vertex_index], 
                verts[indices[3*i+1].vertex_index], verts[indices[3*i+2].vertex_index]));
            SDFace *f = fs.back();
            f->index = fs.size() - 1;

            for (int j = 0; j < 4; j++) // allocate children
                f->child[j] = new SDFace(nullptr, nullptr, nullptr);

            for (int j = 0; j < 3; j++)
            {
                verts[indices[3*i+j].vertex_index]->startFace = f;

                SDEdge e(verts[indices[3*i+j].vertex_index],
                    verts[indices[3*i+NEXT(j)].vertex_index]);
                if (auto prev = edges.find(e); prev == edges.end())
                {
                    e.verts[0]->valence++;
                    e.verts[1]->valence++;
                    e.fs[0] = f;
                    e.f0edgeNum = j;
                    e.index = 1;
                    edges.emplace(e);
                } else
                {
                    // prev->f1edgeNum = j;
                    prev->fs[0]->neighborFs[prev->f0edgeNum] = f;
                    f->neighborFs[j] = (*prev).fs[0];
                    SDEdge newE = *prev;
                    newE.fs[1] = f;
                    newE.f1edgeNum = j;
                    sdshapes[shapei].es.insert(newE);
                    newE.index = prev->index + 1;
                    if (newE.index == 3)
                    {
                        std::cout << "ERROR: model not Valid\n";
                        exit(0);
                    }
                    edges.erase(prev);
                    edges.insert(newE);
                }
            }
        }
        for (SDEdge e: edges)
        {
            sdshapes[shapei].es.insert(e);
        }
        if (shapes[shapei].mesh.smoothing_group_ids[0])
            sdshapes[shapei].smooth = 0;
    }

    std::cout << "model valid\n";
    sdattrib.vs.swap(verts);
    sdattrib.ts.swap(vts);
    sdattrib.ns.swap(vns);
}

/*
* sdattrib: input
* sdshapes: input
* attrib: output
* shapes: output
*/
static void createObjShapes(const SDAttrib &sdattrib, const std::vector<SDShape> &sdshapes, 
    tinyobj::attrib_t &attrib, std::vector<tinyobj::shape_t> &shapes)
{
    // convert to output format
    std::vector<real_t> newVertices;
    std::vector<real_t> newNormals;
    for (std::vector<real_t>::size_type i = 0; i < sdattrib.vs.size(); i++)
    {
        for (int j = 0; j < 3; j++)
            newVertices.emplace_back(sdattrib.vs[i]->p[j]);
    }
    for (std::vector<real_t>::size_type i = 0; i < sdattrib.ns.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            newNormals.emplace_back(sdattrib.ns[i]->p[j]);
        }
    }
    attrib.vertices.swap(newVertices);
    attrib.normals.swap(newNormals);
    // set other attribs to empty
    attrib.texcoords = std::vector<real_t>{};

    for (r_size_type i = 0; i < sdshapes.size(); i++)
    {
        std::vector<tinyobj::index_t> newIndices;
        std::vector<unsigned int> newNumFaceVertices;
        std::vector<unsigned int> newSmoothGroupIds;
        std::vector<int> newMaterialIds;
        for (const SDFace *sdface: sdshapes[i].fs)
        {
            for (int j = 0; j < 3; j++)
            {
                tinyobj::index_t index;
                index.vertex_index = sdface->verts[j]->index;
                index.normal_index = sdface->vns[j]->index;
                index.texcoord_index = -1;
                newIndices.emplace_back(index);
            }
            newNumFaceVertices.emplace_back(3);
            newSmoothGroupIds.emplace_back(0);
            newMaterialIds.emplace_back(-1);
        }
        shapes[i].mesh.indices.swap(newIndices);
        shapes[i].mesh.num_face_vertices.swap(newNumFaceVertices);
        shapes[i].mesh.material_ids.swap(newMaterialIds);
        shapes[i].mesh.smoothing_group_ids.swap(newSmoothGroupIds);
    }
}

static void retainNormalVector(SDAttrib &attrib, SDShape &shape, bool smooth = false)
{
    std::vector<SDNormal*> ns;
    for (SDFace *face: shape.fs)
    {
        // compute normal vector for face
        SDVertex *v0 = face->verts[0], *v1 = face->verts[1], *v2 = face->verts[2];
        real_t a1 = v1->p[0] - v0->p[0];
        real_t a2 = v1->p[1] - v0->p[1];
        real_t a3 = v1->p[2] - v0->p[2];
        
        real_t b1 = v2->p[0] - v1->p[0];
        real_t b2 = v2->p[1] - v1->p[1];
        real_t b3 = v2->p[2] - v1->p[2];

        SDNormal *normal = new SDNormal(a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1);
        real_t len = sqrt(normal->p *normal->p);
        normal->p = normal->p / len;

        ns.push_back(normal);
        normal->index = ns.size() - 1;
        face->vns[0] = face->vns[1] = face->vns[2] = normal;
    }
    attrib.ns.insert(attrib.ns.end(), ns.begin(), ns.end());
}

static void retainNormalVector(SDAttrib &attrib, std::vector<SDShape> &shapes, bool smooth = false)
{
    std::vector<SDNormal*> ns;
    for (SDShape &shape: shapes)
    {
        for (SDFace *face: shape.fs)
        {
            // compute normal vector for _face
            SDVertex *v0 = face->verts[0], *v1 = face->verts[1], *v2 = face->verts[2];
            real_t a1 = v1->p[0] - v0->p[0];
            real_t a2 = v1->p[1] - v0->p[1];
            real_t a3 = v1->p[2] - v0->p[2];
            
            real_t b1 = v2->p[0] - v0->p[0];
            real_t b2 = v2->p[1] - v0->p[1];
            real_t b3 = v2->p[2] - v0->p[2];

            SDNormal *normal = new SDNormal(a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1);
            real_t len = sqrt(normal->p *normal->p);
            normal->p = normal->p / len;

            ns.push_back(normal);
            normal->index = ns.size() - 1;
            face->vns[0] = face->vns[1] = face->vns[2] = normal;
        }
    }
    attrib.ns.swap(ns);
}

/*
* input: v
* output: range begining at d_first
* describe:
*   find all vertices surrounding the vertex v. after process, 
*   vertices will be output to range begining at d_first. 
*   d_first should meet requirements of LegacyOutputIterator
*/
template <class OutputIt>
inline bool findAllNeighbors(const SDVertex *v, OutputIt d_first)
{// get all neighbors of SDVertex v. Return true if v is boundary, false either
    SDFace *f = v->startFace;
    do {
        *d_first = f->verts[NEXT(f->getVertIndex(v))];
        f = f->nextFace(v);
    } while (f && f != v->startFace);

    if (!f) 
    {
        f = v->startFace;
        do {
            *d_first = f->verts[PREV(f->getVertIndex(v))];
            f = f->prevFace(v);
        } while (f && f != v->startFace);
        return true;
    }
    return false;
}


}

#endif