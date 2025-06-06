#ifndef SOLVER_BASE_H
#define SOLVER_BASE_H

#include "tiny_obj_loader.h"
#include <mutex>

struct IntParam
{
    const char *name;
    int value;
};

struct FltParam
{
    const char *name;
    float value;
};

class SolverBase
{
    std::vector<IntParam> intParams;
    std::vector<FltParam> fltParams;
public:
    std::mutex finMtx;
    std::mutex updateMtx;
    std::mutex sourceMtx;
    bool hasFinished = false;
    bool update = false;

    virtual const char *getName() = 0;
    virtual void run(tinyobj::attrib_t&, std::vector<tinyobj::shape_t>&) = 0;
    std::vector<IntParam> getIntParamList() const {return intParams;}
    std::vector<FltParam> getFltParamList() const {return fltParams;}

    void addIntParam(const char *name, int val) {intParams.push_back(IntParam{name, val});}
    void setIntParam(int index, int val) {intParams[index].value = val;}
    int getIntParam(int index) {return intParams[index].value;}

    void addFltParam(const char *name, float val) {fltParams.push_back(FltParam{name, val});}
    void setFltParam(int index, float val) {fltParams[index].value = val;}
    float getFltParam(int index) {return fltParams[index].value;}
};

#endif