// 4D-only grid lookup and incremental output. Public momenta are in Angstrom^-1.
// No expanded LCAO HDF5 files: state indices address the resident first-BZ arrays.
#pragma once
#include <H5Cpp.h>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace spadexp4d {
using Clock = std::chrono::steady_clock;
inline double seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
inline int wrap(long long x, int n) { return static_cast<int>((x % n + n) % n); }
inline long long nearest(double x) { return static_cast<long long>(std::nearbyint(x)); } // numpy rint tie rule
inline double scalar(H5::Group& g, const char* key) { double v; g.openAttribute(key).read(H5::PredType::NATIVE_DOUBLE,&v); return v; }
inline void write_scalar(H5::H5Object& g, const char* key, double v) {
    H5::DataSpace space(H5S_SCALAR); auto a=g.createAttribute(key,H5::PredType::NATIVE_DOUBLE,space);a.write(H5::PredType::NATIVE_DOUBLE,&v);
}
template<class T> inline void dataset(H5::H5File& f,const char* name,const H5::PredType& type,
                                    const std::vector<hsize_t>& shape,const T* data) {
    H5::DataSpace space(static_cast<int>(shape.size()),shape.data());
    auto ds=f.createDataSet(name,type,space);ds.write(data,type);
}
struct BatchGrid {
    int n[3], nx, ny, state_count;
    bool raw_spectrum=false;
    double energy_step=0;
    std::array<std::array<double,3>,3> g{};
    std::vector<double> eb, distance, weights;
    std::vector<int> index;
    double kmin, kmax, h2m, band_eh, hv, v0, de, max_distance=0;
    std::string output;
    Clock::time_point started=Clock::now();
    double kernel_seconds=0;
    BatchGrid(H5::H5File& input, const char* output_file, int xcount, int ycount,
              double photon, double inner, double width):nx(xcount),ny(ycount),hv(photon),v0(inner),de(width),output(output_file) {
        auto mesh=input.openGroup("Mesh");mesh.openAttribute("Divisions").read(H5::PredType::NATIVE_INT,n);
        mesh.openAttribute("ReciprocalVectors_Ainv").read(H5::PredType::NATIVE_DOUBLE,g.data());
        long long count=1;
        for(int i=0;i<3;++i) {if(n[i]<1 || n[i]>256)throw std::runtime_error("Invalid mesh divisions");count*=n[i];}
        state_count=static_cast<int>(count);
        auto group=input.openGroup("Batch");kmin=scalar(group,"KMin");kmax=scalar(group,"KMax");h2m=scalar(group,"H2M");band_eh=scalar(group,"BandHartree_eV");
        raw_spectrum=group.attrExists("RawSpectrum") && scalar(group,"RawSpectrum")!=0;
        if(raw_spectrum){energy_step=scalar(group,"EnergyStep_eV");if(energy_step<=0)throw std::runtime_error("Invalid energy bin width");}
        auto energy=group.openDataSet("BindingEnergies");auto space=energy.getSpace();hsize_t size;
        if(space.getSimpleExtentNdims()!=1)throw std::runtime_error("Invalid energy array");space.getSimpleExtentDims(&size);
        if(size<1 || size>1000000 || nx<2 || ny<2 || kmax<=kmin || h2m<=0)throw std::runtime_error("Invalid batch geometry");
        eb.resize(size);energy.read(eb.data(),H5::PredType::NATIVE_DOUBLE);
        if(std::abs(g[1][0])+std::abs(g[0][2])+std::abs(g[1][2])+std::abs(g[2][0])+std::abs(g[2][1])>1e-12)
            throw std::runtime_error("This nearest-grid search requires triangular in-plane G and normal G3");
        for(double e:eb)if(!std::isfinite(e)||hv+v0-e<=0)throw std::runtime_error("Invalid binding energy");
        index.resize(nx*ny);distance.resize(nx*ny);weights.resize(nx*ny);
    }
    // Exact Cartesian nearest neighbour on the infinite triangular grid.
    // The loop bound follows from the best candidate's Cartesian distance.
    int lookup(double x,double y,double z,double& error)const {
        double sx=g[0][0]/n[0],sy=g[0][1]/n[0],ty=g[1][1]/n[1],sz=g[2][2]/n[2];
        long long cx=nearest(x/sx),iz=nearest(z/sz),bx=cx,by=nearest((y-cx*sy)/ty);
        auto norm=[&](long long ix,long long iy){double dx=ix*sx-x,dy=ix*sy+iy*ty-y,dz=iz*sz-z;return dx*dx+dy*dy+dz*dz;};
        double best=norm(bx,by);int radius=static_cast<int>(std::ceil(std::sqrt(best)/std::abs(sx)));
        for(int shift=-radius;shift<=radius;++shift){long long ix=cx+shift,iy=nearest((y-ix*sy)/ty);double d=norm(ix,iy);
            if(d<best-1e-14){best=d;bx=ix;by=iy;}}
        error=std::sqrt(best);return (wrap(bx,n[0])*n[1]+wrap(by,n[1]))*n[2]+wrap(iz,n[2]);
    }
    void prepare(size_t slice,double (*points)[3],double** band,int nb,double ef,double eh,double au_ang) {
        max_distance=0;double kt2=(hv+v0-eb[slice])/h2m;
        for(int iy=0;iy<ny;++iy)for(int ix=0;ix<nx;++ix){int i=iy*nx+ix;
            double x=kmin+(kmax-kmin)*ix/(nx-1),y=kmin+(kmax-kmin)*iy/(ny-1),z2=kt2-x*x-y*y;
            double z=std::sqrt(std::max(z2,0.0));points[i][0]=x*au_ang;points[i][1]=y*au_ang;points[i][2]=z*au_ang;
            index[i]=lookup(x,y,z,distance[i]);max_distance=std::max(max_distance,distance[i]);weights[i]=0;
            if(z2>=0)for(int b=0;b<nb;++b){
                if(raw_spectrum){
                    double pixel=((band[index[i]][b]-ef)*eh+eb[slice])/energy_step;
                    if(std::floor(pixel+0.5)==0)weights[i]+=1.0/energy_step;
                }else{
                    double delta=((band[index[i]][b]-ef)*band_eh+eb[slice])/de;
                    weights[i]+=std::exp(-0.5*delta*delta);
                }
            }
        }
    }
    std::string slice_path(size_t slice)const {std::ostringstream s;s<<output<<'.'<<std::setw(6)<<std::setfill('0')<<slice<<".h5";return s.str();}
    void save(size_t slice,double** pad,double kernel_time) {
        auto path=slice_path(slice),temporary=path+".tmp";
        {
            H5::H5File f(temporary,H5F_ACC_TRUNC);
            dataset(f,"PAD",H5::PredType::NATIVE_DOUBLE,{static_cast<hsize_t>(ny),static_cast<hsize_t>(nx)},pad[0]);
            dataset(f,"BandWeight",H5::PredType::NATIVE_DOUBLE,{static_cast<hsize_t>(ny),static_cast<hsize_t>(nx)},weights.data());
            dataset(f,"MeshIndex",H5::PredType::NATIVE_INT,{static_cast<hsize_t>(ny),static_cast<hsize_t>(nx)},index.data());
            write_scalar(f,"RawSpectrum",raw_spectrum);
            write_scalar(f,"EB",eb[slice]);write_scalar(f,"MaxSamplingDistance_Ainv",max_distance);
            write_scalar(f,"KernelSeconds",kernel_time);write_scalar(f,"Complete",1);
        }
        std::filesystem::rename(temporary,path);kernel_seconds+=kernel_time;
        std::cout<<"[4D_SLICE] "<<slice<<" EB="<<eb[slice]<<" kernel="<<kernel_time<<" s"<<std::endl;
    }
    void finish() {
        H5::H5File summary(output,H5F_ACC_TRUNC);
        write_scalar(summary,"Complete",1);write_scalar(summary,"SliceCount",eb.size());
        write_scalar(summary,"StateCount",state_count);write_scalar(summary,"CoefficientLoadCount",1);
        write_scalar(summary,"KernelSeconds",kernel_seconds);write_scalar(summary,"BatchWallSeconds",seconds(started));
        std::cout<<"[4D_BATCH] "<<eb.size()<<" slices, coefficients loaded once, wall="<<seconds(started)<<" s"<<std::endl;
    }
};
} // namespace spadexp4d
