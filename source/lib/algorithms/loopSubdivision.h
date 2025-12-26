#ifndef LOOP_SUBDIVISION_H
#define LOOP_SUBDIVISION_H

#include <iterator>
#include <set>
#include <vector>
#include "MeshTopology.h"
#include "solverBase.h"

#ifndef NDEBUG
#include <iostream>
#endif


namespace LoopSubdivision 
{ // ref: https://pbr-book.org/3ed-2018/Shapes/Subdivision_Surfaces#SDVertex::regular
using namespace meshTopology;

static void loopSubdivision(tinyobj::attrib_t&, std::vector<tinyobj::shape_t>&);
inline real_t beta(int valence);
template <class V>
void computeWeightedPosition(const std::vector<V*> &comp, real_t beta, Point &resp);

class LoopSolver: public SolverBase
{
    const char *name = "Loop Subdivision";
public:
    const char *getName() {return name;}
    void run(tinyobj::attrib_t &attrib, std::vector<tinyobj::shape_t> &shapes) 
    {
        std::vector<SDShape> sdshapes;
        SDAttrib sdattrib;
        createSDShapes(attrib, shapes, sdattrib, sdshapes);

        loopSubdivision(sdattrib, sdshapes);

        finMtx.lock();
        hasFinished = true;
        sourceMtx.lock();
        retainNormalVector(sdattrib, sdshapes);
        createObjShapes(sdattrib, sdshapes, attrib, shapes);
        sourceMtx.unlock();
        finMtx.unlock();
    }

    /*
    * Perform Loop subdivision. Save the mesh after subdivision in sdattrib and sdshapes.
        normal vector and textcoord are not neeeded. In output, we just leave sdattrib.ns and sdattrib.ts empty.
    * input & output: sdattrib. sdattrib saves vertex, normal vector and textcoord of the mesh.
        sdshapes. sdshapes saves topology structure of the mesh.
    * ref: https://pbr-book.org/3ed-2018/Shapes/Subdivision_Surfaces
    */
    void loopSubdivision(SDAttrib &sdattrib, std::vector<SDShape> &sdshapes)
    {
        std::vector<SDVertex*> new_vs;
        
        //原模型中的点位置变换
        for(auto x : sdattrib.vs) {
            auto new_p = new SDVertex();
            std::vector<SDVertex*> neibor;
            findVtxNeighbors(x, std::back_inserter(neibor));
            new_p->index = new_vs.size();
            new_p->valence = neibor.size();
            if(isBoundary(x)) {
                new_p->p = neibor.size() > 1 ? 
                    0.75 * x->p + 0.125 * (neibor.back()->p + neibor[0]->p) : x->p;
            } else {
                int k = neibor.size();
                real_t b = beta(k);
                Point sum(0, 0, 0);
                for(auto x : neibor) sum = sum + x->p;
                new_p->p = (1.0 - k * b) * x->p + b * sum;
            }
            x->child = new_p;
            new_vs.emplace_back(new_p);
        }

        //每条边上新增一个点
        for(auto& shape : sdshapes) {
            for(auto e : shape.es) {
                auto new_p = new SDVertex();
                new_p->valence = 2;
                new_p->index = new_vs.size();
                auto v1 = e.verts[0];
                auto v2 = e.verts[1];
                if(isBoundary(e)) {
                    new_p->p = (v1->p + v2->p) / 2;
                } else {
                    auto v3 = e.fs[0]->otherVertex(&e);
                    auto v4 = e.fs[1]->otherVertex(&e);
                    new_p->p = 0.375 * (v1->p + v2->p) + 0.125 * (v3->p + v4->p);
                }
                new_vs.emplace_back(new_p);

                for(int i = 0; i < 2; ++i) {
                    auto fc = e.fs[i];
                    int id1 = fc->getVertIndex(v1), id2 = fc->getVertIndex(v2);
                    int edId = std::abs(id1 - id2) == 1 ?
                        std::min(id1, id2) : std::max(id1, id2);
                    fc->sdverts[edId] = new_p;
                }
            }
        }

        //重构所有的面
        for(auto& shape : sdshapes) {
            std::vector<SDFace*> new_fs;
            for(auto fc : shape.fs) {
                for(int i = 0; i < 3; ++i) {
                    new_fs.emplace_back(new SDFace(
                        fc->verts[i]->child,
                        fc->sdverts[i],
                        fc->sdverts[(i + 2) % 3]
                    ));
                }
                new_fs.emplace_back(new SDFace(
                    fc->sdverts[0],
                    fc->sdverts[1],
                    fc->sdverts[2]
                ));
                delete fc;
            }
            shape.fs.swap(new_fs);
            shape.es.clear();
        }
        
        for(auto x : sdattrib.vs) delete x;
        sdattrib.vs.swap(new_vs);
    }
};

inline real_t beta(int valence)
{
    if (valence == 3) return 3.f / 16.f;
    return 3.f / (8.f * valence);
}

template <class V>
void computeWeightedPosition(const std::vector<V*> &comp, real_t beta, Point &resp)
{
    for (int i = 0; i < 3; i++)
    {
        resp[i] *= (1 - comp.size() * beta);
        for (V *vertp: comp) 
            resp[i] += vertp->p[i] * beta;
    }
}


} // namespace: loopSubdivision

#endif