#ifndef ROOT3_SUBDIVISION_H
#define ROOT3_SUBDIVISION_H

#include <cmath>
#include <iterator>
#include <vector>
#include "solverBase.h"
#include "MeshTopology.h"

#ifndef NDEBUG
#include <iostream>
#endif

namespace root3Subdivision
{
using namespace meshTopology;

static void root3Subdivision(tinyobj::attrib_t&, std::vector<tinyobj::shape_t>&);
template <class It>
inline void oneRingRule(const It begin, const It end, SDVertex *res, real_t w);
inline real_t alpha(int valence);

class Root3Solver: public SolverBase
{
    const char *name = "Root3 Subdivision";
public:
    const char *getName() {return name;}
    void run(tinyobj::attrib_t &attrib, std::vector<tinyobj::shape_t> &shapes) 
    {
        std::vector<SDShape> sdshapes;
        SDAttrib sdattrib;
        createSDShapes(attrib, shapes, sdattrib, sdshapes);

        root3Subdivision(sdattrib, sdshapes);

        finMtx.lock();
        hasFinished = true;
        sourceMtx.lock();
        retainNormalVector(sdattrib, sdshapes);
        createObjShapes(sdattrib, sdshapes, attrib, shapes);
        sourceMtx.unlock();
        finMtx.unlock();
    }

    /*
    * Perform Root3 subdivision. Save the mesh after subdivision in sdattrib and sdshapes.
        normal vector and textcoord are not needed. In output, we just leave sdattrib.ns and sdattrib.ts empty.
    * input & output: sdattrib. sdattrib saves vertex, normal vector and textcoord of the mesh.
        sdshapes. sdshapes saves topology structure of the mesh.
    * ref: https://dl.acm.org/doi/10.1145/344779.344835
    */
    void root3Subdivision(SDAttrib &sdattrib, std::vector<SDShape> &sdshapes)
    {
        std::vector<SDVertex*> new_vs;
        for(auto x : sdattrib.vs) {
            auto new_p = new_vs.emplace_back(new SDVertex(x->p));
            std::vector<SDVertex*> neibor;
            findVtxNeighbors(x, std::back_inserter(neibor));
            new_p->index = new_vs.size() - 1;
            new_p->valence = neibor.size();
            oneRingRule(neibor.begin(), neibor.end(), new_p, alpha(new_p->valence));
            x->child = new_p;
        }

        for(auto &shape : sdshapes) {
            for(auto fc : shape.fs) {
                auto new_p = new_vs.emplace_back(new SDVertex());
                new_p->index = new_vs.size() - 1;
                new_p->valence = 6;
                for(int i = 0; i < 3; ++i) new_p->p = new_p->p + fc->verts[i]->p;
                new_p->p = new_p->p / 3;
                fc->sdverts[0] = new_p;
            }
        }

        for(auto &shape : sdshapes) {
            std::vector<SDFace*> new_fs;
            for(auto fc : shape.fs) {
                for(int i = 0; i < 3; ++i) {
                    auto nfc = fc->neighborFs[i];
                    if(nfc) {
                        auto v = fc->sdverts[0];
                        auto nv = nfc->sdverts[0];
                        auto p = fc->verts[i]->child;
                        new_fs.emplace_back(new SDFace(
                            fc->sdverts[0],
                            nfc->sdverts[0],
                            fc->verts[i]->child
                        ));
                    }
                }
            }
            for(auto fc : shape.fs) delete fc;
            shape.fs.swap(new_fs);
        }
        for(auto x : sdattrib.vs) delete x;
        for(auto x : sdattrib.ns) delete x;
        for(auto x : sdattrib.ts) delete x;
        sdattrib.vs.swap(new_vs);
        sdattrib.ns.clear();
        sdattrib.ts.clear();
    }
};

template <class It>
inline void oneRingRule(const It begin, const It end, SDVertex *res, real_t w)
{
    real_t wres = 1 - res->valence * w;
    for (int i = 0; i < 3; i++)
        res->p[i] *= wres;
    for (It ele = begin; ele != end; ++ele)
    {
        for (int i: {0, 1, 2})
            res->p[i] += w * (*ele)->p[i];
    }
}

inline real_t alpha(int valence)
{
    switch (valence)
    {
    case 1:
        return 2.f / 9.f;
    case 2:
        return 1.f / 3.f;
    case 3:
        return 5.f / 27.f;
    case 4:
        return 1.f / 9.f;
    case 5:
        break;
    case 6:
        return 1.f / 18.f;    
    default:
        break;
    }
    const float pi = 3.141592653589;
    return (4.f-2.f*std::cos(2.f*pi/valence)) / (9.f*(float)valence);
}

}

#endif