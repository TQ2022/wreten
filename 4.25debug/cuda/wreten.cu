#include <list>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <mechsys/flbm/Domain.h>
#include <mechsys/dem/domain.h>
#include <mechsys/util/util.h>
struct UserData
{

    thrust::device_vector<real> Xmin;
    thrust::device_vector<real> Xmax;
    thrust::device_vector<real> Ymin;
    thrust::device_vector<real> Ymax;
    thrust::device_vector<real> Zmin;
    thrust::device_vector<real> Zmax;

    real *pXmin;
    real *pXmax;
    real *pYmin;
    real *pYmax;
    real *pZmin;
    real *pZmax;

    real ome;
    size_t block;
    real Head;
    real Orig;
    real Tf;
    real dtOut;
    real time;
    real rho;
    Vec3_t Dp;
    std::ofstream oss_ss;
};

// Enumeration to define different boundary condition types.
enum BoundaryConditionType
{
    BCT_XMIN0,
    BCT_XMIN1,
    BCT_XMAX0,
    BCT_XMAX1,
    BCT_YMIN0,
    BCT_YMIN1,
    BCT_YMAX0,
    BCT_YMAX1,
    BCT_ZMIN0,
    BCT_ZMIN1,
    BCT_ZMAX0,
    BCT_ZMAX1
};

__global__ void SetupBoundaryConditions(real *rhoBC, bool *IsSolid, real *F, real3 *Vel, real *Rho, FLBM::lbm_aux *lbmaux, BoundaryConditionType bcType)
{

    int ic = threadIdx.x + blockIdx.x * blockDim.x;
    int Nx = lbmaux[0].Nx;
    int Ny = lbmaux[0].Ny;
    int Nz = lbmaux[0].Nz;

    size_t ib, ix, iy, iz, iv;

    switch (bcType)
    {
    case BCT_XMIN0:

        if (ic >= Ny * Nz)
            return;
        ix = 0;
        iy = ic % Ny;
        iz = (ic / Ny) % Nz;
        ib = ix + iy * Nx + iz * Nx * Ny;
        //if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = &F[iv];
            f[1] = 1.0 / 3.0 * (-2 * f[0] - 4 * f[10] - 4 * f[12] - 4 * f[14] - f[2] - 2 * f[3] - 2 * f[4] - 2 * f[5] - 2 * f[6] - 4 * f[8] + 2 * rhoBC[0]);
            f[7] = 1.0 / 24.0 * (-2 * f[0] - 4 * f[10] - 4 * f[12] - 4 * f[14] - 4 * f[2] + f[3] - 5 * f[4] + f[5] - 5 * f[6] + 20 * f[8] + 2 * rhoBC[0]);
            f[9] = 1.0 / 24.0 * (-2 * f[0] + 20 * f[10] - 4 * f[12] - 4 * f[14] - 4 * f[2] + f[3] - 5 * f[4] - 5 * f[5] + f[6] - 4 * f[8] + 2 * rhoBC[0]);
            f[11] = 1.0 / 24.0 * (-2 * f[0] - 4 * f[10] + 20 * f[12] - 4 * f[14] - 4 * f[2] - 5 * f[3] + f[4] + f[5] - 5 * f[6] - 4 * f[8] + 2 * rhoBC[0]);
            f[13] = 1.0 / 24.0 * (-2 * f[0] - 4 * f[10] - 4 * f[12] + 20 * f[14] - 4 * f[2] - 5 * f[3] + f[4] - 5 * f[5] + f[6] - 4 * f[8] + 2 * rhoBC[0]);

            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {
                //if ((iz==Nz/2)&&(iy==Ny/2)) printf("F %g %g %d \n",F[iv + k],f[k],k);
                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            //  printf("xmin0  rhoBC[0] = %f, Rho[ib] = %f, f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //    rhoBC[0], Rho[ib], f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);


    // printf("xmin0  rhoBC[0] = %f, Rho[ib+Ny*Nz] = %f\n", rhoBC[0], Rho[ib+Ny*Nz]);

            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];

        }

        break;
    case BCT_XMIN1:
        if (ic >= Ny * Nz)
            return;
        ix = 0;
        iy = ic % Ny;
        iz = (ic / Ny) % Nz;
        ib = ix + iy * Nx + iz * Nx * Ny + Nx * Ny * Nz;
        // ib = ix * Nz * Ny + iy * Nz + iz + Nx * Ny * Nz;

        //if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;
            f[1] = 1.0 / 3.0 * (-2 * f[0] - 4 * f[10] - 4 * f[12] - 4 * f[14] - f[2] - 2 * f[3] - 2 * f[4] - 2 * f[5] - 2 * f[6] - 4 * f[8] + 2 * rhoBC[1]);
            f[7] = 1.0 / 24.0 * (-2 * f[0] - 4 * f[10] - 4 * f[12] - 4 * f[14] - 4 * f[2] + f[3] - 5 * f[4] + f[5] - 5 * f[6] + 20 * f[8] + 2 * rhoBC[1]);
            f[9] = 1.0 / 24.0 * (-2 * f[0] + 20 * f[10] - 4 * f[12] - 4 * f[14] - 4 * f[2] + f[3] - 5 * f[4] - 5 * f[5] + f[6] - 4 * f[8] + 2 * rhoBC[1]);
            f[11] = 1.0 / 24.0 * (-2 * f[0] - 4 * f[10] + 20 * f[12] - 4 * f[14] - 4 * f[2] - 5 * f[3] + f[4] + f[5] - 5 * f[6] - 4 * f[8] + 2 * rhoBC[1]);
            f[13] = 1.0 / 24.0 * (-2 * f[0] - 4 * f[10] - 4 * f[12] + 20 * f[14] - 4 * f[2] - 5 * f[3] + f[4] - 5 * f[5] + f[6] - 4 * f[8] + 2 * rhoBC[1]);

            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            // printf("xmin0  rhoBC[1] = %f, Rho[ib] = %f, f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //        rhoBC[1], Rho[ib], f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);
            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];
        }
        break;

    case BCT_XMAX0:
        if (ic >= Ny * Nz)
            return;
        ix = Nx - 1;
        iy = ic % Ny;
        iz = ic / Ny;
        // idx = ix + iy * Nx + iz * Nx * Ny;

        ib = ix + iy * Nx + iz * Nx * Ny;

        //if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;
            f[2] = 1 / 3.0 * (-2 * f[0] - f[1] - 2 * (2 * f[11] + 2 * f[13] + f[3] + f[4] + f[5] + f[6] + 2 * f[7] + 2 * f[9] - rhoBC[0]));
            f[8] = 1 / 24.0 * (-2 * f[0] - 4 * f[1] - 4 * f[11] - 4 * f[13] - 5 * f[3] + f[4] - 5 * f[5] + f[6] + 20 * f[7] - 4 * f[9] + 2 * rhoBC[0]);
            f[10] = 1 / 24.0 * (-2 * f[0] - 4 * f[1] - 4 * f[11] - 4 * f[13] - 5 * f[3] + f[4] + f[5] - 5 * f[6] - 4 * f[7] + 20 * f[9] + 2 * rhoBC[0]);
            f[12] = 1 / 24.0 * (-2 * f[0] - 4 * f[1] + 20 * f[11] - 4 * f[13] + f[3] - 5 * f[4] - 5 * f[5] + f[6] - 4 * f[7] - 4 * f[9] + 2 * rhoBC[0]);
            f[14] = 1 / 24.0 * (-2 * f[0] - 4 * f[1] - 4 * f[11] + 20 * f[13] + f[3] - 5 * f[4] + f[5] - 5 * f[6] - 4 * f[7] - 4 * f[9] + 2 * rhoBC[0]);

            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            // printf("rhoBC[0] = %f, Rho[ib] = %f, f[0] = %f, f[2] = %f, f[8] = %f, f[10] = %f, f[12] = %f, f[14] = %f\n", rhoBC[0], Rho[ib], f[0], f[2], f[8], f[10], f[12], f[14]);

            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];

        }

        break;
    case BCT_XMAX1:
        if (ic >= Ny * Nz)
            return;
        ix = Nx - 1;
        iy = ic % Ny;
        iz = ic / Ny;
        // idx = ix + iy * Nx + iz * Nx * Ny + Nx * Ny * Nz;

        ib = ix + iy * Nx + iz * Nx * Ny + Nx * Ny * Nz;

        //if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;
            f[2] = 1 / 3.0 * (-2 * f[0] - f[1] - 2 * (2 * f[11] + 2 * f[13] + f[3] + f[4] + f[5] + f[6] + 2 * f[7] + 2 * f[9] - rhoBC[1]));
            f[8] = 1 / 24.0 * (-2 * f[0] - 4 * f[1] - 4 * f[11] - 4 * f[13] - 5 * f[3] + f[4] - 5 * f[5] + f[6] + 20 * f[7] - 4 * f[9] + 2 * rhoBC[1]);
            f[10] = 1 / 24.0 * (-2 * f[0] - 4 * f[1] - 4 * f[11] - 4 * f[13] - 5 * f[3] + f[4] + f[5] - 5 * f[6] - 4 * f[7] + 20 * f[9] + 2 * rhoBC[1]);
            f[12] = 1 / 24.0 * (-2 * f[0] - 4 * f[1] + 20 * f[11] - 4 * f[13] + f[3] - 5 * f[4] - 5 * f[5] + f[6] - 4 * f[7] - 4 * f[9] + 2 * rhoBC[1]);
            f[14] = 1 / 24.0 * (-2 * f[0] - 4 * f[1] - 4 * f[11] + 20 * f[13] + f[3] - 5 * f[4] + f[5] - 5 * f[6] - 4 * f[7] - 4 * f[9] + 2 * rhoBC[1]);

            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }

            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];
        }
        break;
    case BCT_YMIN0:
        if (ic >= Nx * Nz)
            return;
        ix = ic % Nx;
        iy = 0;
        iz = ic / Nx;
        // idx = ix + iy * Nx + iz * Nx * Ny;
        ib = ix + iy * Nx + iz * Nx * Ny;

        if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;

            f[3] = 1 / 3.0 * (-2 * f[0] - 2 * f[1] - 4 * f[10] - 4 * f[11] - 4 * f[13] - 2 * f[2] - f[4] - 2 * f[5] - 2 * f[6] - 4 * f[8] + 2 * rhoBC[0]);
            f[7] = 1 / 24.0 * (-2 * f[0] + f[1] - 4 * f[10] - 4 * f[11] - 4 * f[13] - 5 * f[2] - 4 * f[4] + f[5] - 5 * f[6] + 20 * f[8] + 2 * rhoBC[0]);
            f[9] = 1 / 24.0 * (-2 * f[0] + f[1] + 20 * f[10] - 4 * f[11] - 4 * f[13] - 5 * f[2] - 4 * f[4] - 5 * f[5] + f[6] - 4 * f[8] + 2 * rhoBC[0]);
            f[12] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[10] + 20 * f[11] - 4 * f[13] + f[2] - 4 * f[4] - 5 * f[5] + f[6] - 4 * f[8] + 2 * rhoBC[0]);
            f[14] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[10] - 4 * f[11] + 20 * f[13] + f[2] - 4 * f[4] + f[5] - 5 * f[6] - 4 * f[8] + 2 * rhoBC[0]);

            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            // printf("YMIN0  rhoBC[0] = %f, Rho[ib] = %f, f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //        rhoBC[0], Rho[ib], f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);
            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];
        }
        break;
    case BCT_YMIN1:
        if (ic >= Nx * Nz)
            return;
        ix = ic % Nx;
        iy = 0;
        iz = ic / Nx;
        ib = ix + iy * Nx + iz * Nx * Ny + Nx * Ny * Nz;

        if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;
            f[3] = 1 / 3.0 * (-2 * f[0] - 2 * f[1] - 4 * f[10] - 4 * f[11] - 4 * f[13] - 2 * f[2] - f[4] - 2 * f[5] - 2 * f[6] - 4 * f[8] + 2 * rhoBC[1]);
            f[7] = 1 / 24.0 * (-2 * f[0] + f[1] - 4 * f[10] - 4 * f[11] - 4 * f[13] - 5 * f[2] - 4 * f[4] + f[5] - 5 * f[6] + 20 * f[8] + 2 * rhoBC[1]);
            f[9] = 1 / 24.0 * (-2 * f[0] + f[1] + 20 * f[10] - 4 * f[11] - 4 * f[13] - 5 * f[2] - 4 * f[4] - 5 * f[5] + f[6] - 4 * f[8] + 2 * rhoBC[1]);
            f[12] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[10] + 20 * f[11] - 4 * f[13] + f[2] - 4 * f[4] - 5 * f[5] + f[6] - 4 * f[8] + 2 * rhoBC[1]);
            f[14] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[10] - 4 * f[11] + 20 * f[13] + f[2] - 4 * f[4] + f[5] - 5 * f[6] - 4 * f[8] + 2 * rhoBC[1]);
            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            // printf("YMIN1  rhoBC[1] = %f, Rho[ib] = %f, F + iv,f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //        rhoBC[1], Rho[ib], F[iv ],f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);
            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];
        }
        break;
    case BCT_YMAX0:
        if (ic >= Nx * Nz)
            return;
        ix = ic % Nx;
        iy = Ny - 1;
        iz = ic / Nx;

        ib = ix + iy * Nx + iz * Nx * Ny;

        if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;
            f[4] = 1 / 3.0 * (-2 * f[0] - 2 * f[1] - 4 * f[12] - 4 * f[14] - 2 * f[2] - f[3] - 2 * f[5] - 2 * f[6] - 4 * f[7] - 4 * f[9] + 2 * rhoBC[0]);
            f[8] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[12] - 4 * f[14] + f[2] - 4 * f[3] - 5 * f[5] + f[6] + 20 * f[7] - 4 * f[9] + 2 * rhoBC[0]);
            f[10] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[12] - 4 * f[14] + f[2] - 4 * f[3] + f[5] - 5 * f[6] - 4 * f[7] + 20 * f[9] + 2 * rhoBC[0]);
            f[11] = 1 / 24.0 * (-2 * f[0] + f[1] + 20 * f[12] - 4 * f[14] - 5 * f[2] - 4 * f[3] + f[5] - 5 * f[6] - 4 * f[7] - 4 * f[9] + 2 * rhoBC[0]);
            f[13] = 1 / 24.0 * (-2 * f[0] + f[1] - 4 * f[12] + 20 * f[14] - 5 * f[2] - 4 * f[3] - 5 * f[5] + f[6] - 4 * f[7] - 4 * f[9] + 2 * rhoBC[0]);
            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            // printf("YMMAX0  rhoBC[0] = %f, Rho[ib] = %f, f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //        rhoBC[0], Rho[ib], f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);
            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];
        }
        break;
    case BCT_YMAX1:
        if (ic >= Nx * Nz)
            return;
        ix = ic % Nx;
        iy = Ny - 1;
        iz = ic / Nx;
        // idx = ix + iy * Nx + iz * Nx * Ny + Nx * Ny * Nz;
        ib = ix + iy * Nx + iz * Nx * Ny + lbmaux[0].Nx * lbmaux[0].Ny * lbmaux[0].Nz;

        if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;
            f[4] = 1 / 3.0 * (-2 * f[0] - 2 * f[1] - 4 * f[12] - 4 * f[14] - 2 * f[2] - f[3] - 2 * f[5] - 2 * f[6] - 4 * f[7] - 4 * f[9] + 2 * rhoBC[1]);
            f[8] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[12] - 4 * f[14] + f[2] - 4 * f[3] - 5 * f[5] + f[6] + 20 * f[7] - 4 * f[9] + 2 * rhoBC[1]);
            f[10] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[12] - 4 * f[14] + f[2] - 4 * f[3] + f[5] - 5 * f[6] - 4 * f[7] + 20 * f[9] + 2 * rhoBC[1]);
            f[11] = 1 / 24.0 * (-2 * f[0] + f[1] + 20 * f[12] - 4 * f[14] - 5 * f[2] - 4 * f[3] + f[5] - 5 * f[6] - 4 * f[7] - 4 * f[9] + 2 * rhoBC[1]);
            f[13] = 1 / 24.0 * (-2 * f[0] + f[1] - 4 * f[12] + 20 * f[14] - 5 * f[2] - 4 * f[3] - 5 * f[5] + f[6] - 4 * f[7] - 4 * f[9] + 2 * rhoBC[1]);
            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            // printf("Ymax1  rhoBC[1] = %f, Rho[ib] = %f, f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //        rhoBC[1], Rho[ib], f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);
            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];
        }
        break;
    case BCT_ZMIN0:
        if (ic >= Nx * Ny)
            return;
        ix = ic % Nx;
        iy = (ic / Nx) % Ny;
        iz = 0;
        //  idx = ix + iy * Nx + iz * Nx * Ny;
        ib = ix + iy * Nx + iz * Nx * Ny;

        if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;
            f[5] = 1 / 3.0 * (-2 * f[0] - 2 * f[1] - 4 * f[12] - 4 * f[13] - 2 * f[2] - 2 * f[3] - 2 * f[4] - f[6] - 4 * f[8] - 4 * f[9] + 2 * rhoBC[0]);
            f[7] = 1 / 24.0 * (-2 * f[0] + f[1] - 4 * f[12] - 4 * f[13] - 5 * f[2] + f[3] - 5 * f[4] - 4 * f[6] + 20 * f[8] - 4 * f[9] + 2 * rhoBC[0]);
            f[10] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[12] - 4 * f[13] + f[2] - 5 * f[3] + f[4] - 4 * f[6] - 4 * f[8] + 20 * f[9] + 2 * rhoBC[0]);
            f[11] = 1 / 24.0 * (-2 * f[0] + f[1] + 20 * f[12] - 4 * f[13] - 5 * f[2] - 5 * f[3] + f[4] - 4 * f[6] - 4 * f[8] - 4 * f[9] + 2 * rhoBC[0]);
            f[14] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[12] + 20 * f[13] + f[2] + f[3] - 5 * f[4] - 4 * f[6] - 4 * f[8] - 4 * f[9] + 2 * rhoBC[0]);

            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            // printf("ZMIN)  rhoBC[0] = %f, Rho[ib] = %f, f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //        rhoBC[0], Rho[ib], f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);
            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];
        }

        break;
    case BCT_ZMIN1:
        if (ic >= Nx * Ny)
            return;
        ix = ic % Nx;
        iy = (ic / Nx) % Ny;
        iz = 0;

        ib = ix + iy * Nx + iz * Nx * Ny + Nx * Ny * Nz;

        if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;
            f[5] = 1 / 3.0 * (-2 * f[0] - 2 * f[1] - 4 * f[12] - 4 * f[13] - 2 * f[2] - 2 * f[3] - 2 * f[4] - f[6] - 4 * f[8] - 4 * f[9] + 2 * rhoBC[1]);
            f[7] = 1 / 24.0 * (-2 * f[0] + f[1] - 4 * f[12] - 4 * f[13] - 5 * f[2] + f[3] - 5 * f[4] - 4 * f[6] + 20 * f[8] - 4 * f[9] + 2 * rhoBC[1]);
            f[10] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[12] - 4 * f[13] + f[2] - 5 * f[3] + f[4] - 4 * f[6] - 4 * f[8] + 20 * f[9] + 2 * rhoBC[1]);
            f[11] = 1 / 24.0 * (-2 * f[0] + f[1] + 20 * f[12] - 4 * f[13] - 5 * f[2] - 5 * f[3] + f[4] - 4 * f[6] - 4 * f[8] - 4 * f[9] + 2 * rhoBC[1]);
            f[14] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[12] + 20 * f[13] + f[2] + f[3] - 5 * f[4] - 4 * f[6] - 4 * f[8] - 4 * f[9] + 2 * rhoBC[1]);

            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            // printf("ZMIN1  rhoBC[1 = %f, Rho[ib] = %f, f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //        rhoBC[0], Rho[ib], f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);
            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];

            // printf("zmin1   f[0]= %f , f[1] =%f, f[2] =%f, f[3] =%f, f[4] =%f, f[5] =%f, f[6] =%f, f[7] =%f, f[8] =%f, f[9] =%f, f[10] =%f, f[11] =%f, f[12] =%f, f[13] =%f, f[14] =%f\n", f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);
        }
        break;
    case BCT_ZMAX0:
        if (ic >= Nx * Ny)
            return;
        ix = ic % Nx;
        iy = (ic / Nx) % Ny;
        iz = Nz - 1;

        ib = ix + iy * Nx + iz * Nx * Ny;

        if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;
            f[6] = 1 / 3.0 * (-2 * f[0] - 2 * f[1] - 4 * f[10] - 4 * f[11] - 4 * f[14] - 2 * f[2] - 2 * f[3] - 2 * f[4] - f[5] - 4 * f[7] + 2 * rhoBC[0]);
            f[8] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[10] - 4 * f[11] - 4 * f[14] + f[2] - 5 * f[3] + f[4] - 4 * f[5] + 20 * f[7] + 2 * rhoBC[0]);
            f[9] = 1 / 24.0 * (-2 * f[0] + f[1] + 20 * f[10] - 4 * f[11] - 4 * f[14] - 5 * f[2] + f[3] - 5 * f[4] - 4 * f[5] - 4 * f[7] + 2 * rhoBC[0]);
            f[12] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[10] + 20 * f[11] - 4 * f[14] + f[2] + f[3] - 5 * f[4] - 4 * f[5] - 4 * f[7] + 2 * rhoBC[0]);
            f[13] = 1 / 24.0 * (-2 * f[0] + f[1] - 4 * f[10] - 4 * f[11] + 20 * f[14] - 5 * f[2] - 5 * f[3] + f[4] - 4 * f[5] - 4 * f[7] + 2 * rhoBC[0]);
            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            //  printf("ZMAX0  rhoBC[0] = %f, Rho[ib] = %f, F + iv,f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //                    rhoBC[0], Rho[ib], F[iv ],f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);

            // printf("ZMax0  rhoBC[0] = %f, Rho[ib] = %f, f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //        rhoBC[0], Rho[ib], f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);
            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];
        }
        break;

    case BCT_ZMAX1:
        if (ic >= Nx * Ny)
            return;
        ix = ic % Nx;
        iy = (ic / Nx) % Ny;
        iz = Nz - 1;
        ib = ix + iy * Nx + iz * Nx * Ny + Nx * Ny * Nz;
        // idx = ix + iy * Nx + iz * Nx * Ny + Nx * Ny * Nz; // this is the idx for fluid 1 and there is no lbmaux[1]

        if (!IsSolid[ib])
        {
            size_t iv = ib * lbmaux[0].Nneigh;
            real *f = F + iv;
            f[6] = 1 / 3.0 * (-2 * f[0] - 2 * f[1] - 4 * f[10] - 4 * f[11] - 4 * f[14] - 2 * f[2] - 2 * f[3] - 2 * f[4] - f[5] - 4 * f[7] + 2 * rhoBC[1]);
            f[8] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[10] - 4 * f[11] - 4 * f[14] + f[2] - 5 * f[3] + f[4] - 4 * f[5] + 20 * f[7] + 2 * rhoBC[1]);
            f[9] = 1 / 24.0 * (-2 * f[0] + f[1] + 20 * f[10] - 4 * f[11] - 4 * f[14] - 5 * f[2] + f[3] - 5 * f[4] - 4 * f[5] - 4 * f[7] + 2 * rhoBC[1]);
            f[12] = 1 / 24.0 * (-2 * f[0] - 5 * f[1] - 4 * f[10] + 20 * f[11] - 4 * f[14] + f[2] + f[3] - 5 * f[4] - 4 * f[5] - 4 * f[7] + 2 * rhoBC[1]);
            f[13] = 1 / 24.0 * (-2 * f[0] + f[1] - 4 * f[10] - 4 * f[11] + 20 * f[14] - 5 * f[2] - 5 * f[3] + f[4] - 4 * f[5] - 4 * f[7] + 2 * rhoBC[1]);
            Rho[ib] = 0.0;
            Vel[ib] = make_real3(0.0, 0.0, 0.0);
            for (size_t k = 0; k < lbmaux[0].Nneigh; k++)
            {

                Rho[ib] += F[iv + k];
                Vel[ib] = Vel[ib] + F[iv + k] * lbmaux[0].C[k];
            }
            // printf("ZMAX1  rhoBC[1] = %f, Rho[ib] = %f, f[0] = %f, f[1] = %f, f[2] = %f, f[3] = %f, f[4] = %f, f[5] = %f, f[6] = %f, f[7] = %f, f[8] = %f, f[9] = %f, f[10] = %f, f[11] = %f, f[12] = %f, f[13] = %f, f[14] = %f\n",
            //        rhoBC[1], Rho[ib], f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14]);
            Vel[ib] = (lbmaux[0].Cs / Rho[ib]) * Vel[ib];
        }
        break;
    }
}
void Setup(FLBM::Domain &dom, void *UD)
{

    UserData &dat = (*static_cast<UserData *>(UD));
    if (dom.Time > dat.time)
    {
        dat.time += dat.dtOut;
    }
    real a = M_PI / dat.ome;
    real rho = dat.Head * ((1.0 / a) * (dat.time - a * (floor(dat.time / a) + 0.5)) * pow(-1.0, floor(dat.time / a)) + 0.5) + dat.Orig;
    real rho0min;
    real rho1min;
    real rho0max;
    real rho1max;
    int sizeX = dom.Ndim(0);
    int sizeY = dom.Ndim(1);
    int sizeZ = dom.Ndim(2);

    real *pXmin = thrust::raw_pointer_cast(dat.Xmin.data());
    real *pXmax = thrust::raw_pointer_cast(dat.Xmax.data());
    real *pYmin = thrust::raw_pointer_cast(dat.Ymin.data());
    real *pYmax = thrust::raw_pointer_cast(dat.Ymax.data());
    real *pZmin = thrust::raw_pointer_cast(dat.Zmin.data());
    real *pZmax = thrust::raw_pointer_cast(dat.Zmax.data());

    if (fabs(dat.Dp(0)) > 1.0e-12)
    {
        rho0min = 0.999 * ((rho - dat.rho) * dat.Dp(0) + dat.rho);
        rho1min = 0.001 * ((rho - dat.rho) * dat.Dp(0) + dat.rho);
        rho0max = 0.001 * ((dat.rho - rho) * dat.Dp(0) + dat.rho);
        rho1max = 0.999 * ((dat.rho - rho) * dat.Dp(0) + dat.rho);

        // dat.Xmin[0] = 20.0;

        dat.Xmin[0] = rho0min;

        SetupBoundaryConditions<<<(dom.Ndim(1) * dom.Ndim(2)) / dom.Nthread + 1, dom.Nthread>>>(pXmin, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_XMIN0);
        cudaDeviceSynchronize();

        dat.Xmax[0] = rho0max;
        SetupBoundaryConditions<<<(dom.Ndim(1) * dom.Ndim(2)) / dom.Nthread + 1, dom.Nthread>>>(pXmax, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_XMAX0);
        cudaDeviceSynchronize();

        dat.Xmin[1] = rho1min;
        SetupBoundaryConditions<<<(dom.Ndim(1) * dom.Ndim(2)) / dom.Nthread + 1, dom.Nthread>>>(pXmin, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_XMIN1);
        cudaDeviceSynchronize();

        dat.Xmax[1] = rho1max;
        SetupBoundaryConditions<<<(dom.Ndim(1) * dom.Ndim(2)) / dom.Nthread + 1, dom.Nthread>>>(pXmax, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_XMAX1);
        cudaDeviceSynchronize();
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(error));
        }
    }

    if (fabs(dat.Dp(1)) > 1.0e-12)
    {
        dat.Ymin[0] = rho0min = 0.999 * ((rho - dat.rho) * dat.Dp(1) + dat.rho);
        dat.Ymin[1] = rho1min = 0.001 * ((rho - dat.rho) * dat.Dp(1) + dat.rho);
        dat.Ymax[0] = rho0max = 0.001 * ((dat.rho - rho) * dat.Dp(1) + dat.rho);
        dat.Ymax[1] = rho1max = 0.999 * ((dat.rho - rho) * dat.Dp(1) + dat.rho);

        SetupBoundaryConditions<<<dom.Ndim(0) * dom.Ndim(2) / dom.Nthread + 1, dom.Nthread>>>(pYmin, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_YMIN0);
        cudaDeviceSynchronize();

        SetupBoundaryConditions<<<dom.Ndim(0) * dom.Ndim(2) / dom.Nthread + 1, dom.Nthread>>>(pYmax, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_YMAX0);
        cudaDeviceSynchronize();

        SetupBoundaryConditions<<<dom.Ndim(0) * dom.Ndim(2) / dom.Nthread + 1, dom.Nthread>>>(pYmin, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_YMIN1);
        cudaDeviceSynchronize();

        SetupBoundaryConditions<<<dom.Ndim(0) * dom.Ndim(2) / dom.Nthread + 1, dom.Nthread>>>(pYmax, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_YMAX1);
        cudaDeviceSynchronize();
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(error));
        }
    }

    if (fabs(dat.Dp(2)) > 1.0e-12)
    {
        dat.Zmin[0] = rho0min = 0.999 * ((rho - dat.rho) * dat.Dp(2) + dat.rho);
        dat.Zmin[1] = rho1min = 0.001 * ((rho - dat.rho) * dat.Dp(2) + dat.rho);
        dat.Zmax[0] = rho0max = 0.001 * ((dat.rho - rho) * dat.Dp(2) + dat.rho);
        dat.Zmax[1] = rho1max = 0.999 * ((dat.rho - rho) * dat.Dp(2) + dat.rho);

        SetupBoundaryConditions<<<dom.Ndim(0) * dom.Ndim(1) / dom.Nthread + 1, dom.Nthread>>>(pZmin, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_ZMIN0);
        cudaDeviceSynchronize();

        SetupBoundaryConditions<<<dom.Ndim(0) * dom.Ndim(1) / dom.Nthread + 1, dom.Nthread>>>(pZmax, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_ZMAX0);
        cudaDeviceSynchronize();

        SetupBoundaryConditions<<<dom.Ndim(0) * dom.Ndim(1) / dom.Nthread + 1, dom.Nthread>>>(pZmin, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_ZMIN1);
        cudaDeviceSynchronize();

        SetupBoundaryConditions<<<dom.Ndim(0) * dom.Ndim(1) / dom.Nthread + 1, dom.Nthread>>>(pZmax, dom.pIsSolid, dom.pF, dom.pVel, dom.pRho, dom.plbmaux, BCT_ZMAX1);
        cudaDeviceSynchronize();
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(error));
        }
    }
}

void Report(FLBM::Domain &dom, void *UD)
{
    UserData &dat = (*static_cast<UserData *>(UD));
    double water = 0.0;
    double oil = 0.0;
    double Sr = 0.0;
    size_t nw = 0;
    size_t no = 0;

    for (size_t idx = 0; idx < dom.Ncells; idx++)
    {
        iVec3_t coord;
        FLBM::idx2Pt(idx, coord, dom.Ndim);
        bool isSolid = dom.IsSolid[0][coord(0)][coord(1)][coord(2)];
        if (!isSolid)
        {
            double wr = dom.Rho[1][coord(0)][coord(1)][coord(2)];
            double ar = dom.Rho[0][coord(0)][coord(1)][coord(2)];
            if (wr > 0.5 * dat.rho)
            {
                Sr += 1.0;
                water += (wr + ar + dom.Gmix * wr * ar) / 3.0;
                nw++;
            }
            if (ar > 0.5 * dat.rho)
            {
                oil += (wr + ar + dom.Gmix * wr * ar) / 3.0;
                no++;
            }
        }
    }

    double Sf = 0.0;

    for (size_t x = 0; x < dom.Ndim(0); x++)
    {
        for (size_t y = 0; y < dom.Ndim(1); y++)
        {
            for (size_t z = 0; z < dom.Ndim(2); z++)
            {
                if (dom.IsSolid[0][x][y][z])
                {
                    Sf += 1.0;
                }
            }
        }
    } // Solid fraction，即固体分数

    double aaa = Sf / dom.Ncells;

    Sr /= dom.Ncells * (1.0 - aaa);
    if (nw > 0)
        water /= nw;
    if (no > 0)
        oil /= no; //  果然是这里计算有问题

    //   oil = no;
    double rhow = 0.0;
    double rhoo = 0.0;
    size_t nfb = 0;
    size_t nfo = 0;

    for (size_t i = 0; i < dom.Ndim(1); ++i)
    {
        for (size_t j = 0; j < dom.Ndim(2); ++j)
        {
            // 检查边界的固体状态
            if (!dom.IsSolid[0][1][i][j])
            {
                double rho = dom.Rho[0][1][i][j]; //  0 phase  0.999*2
                rhow += rho;
                nfb++;
          //   printf("dom.Rho[0][1][%lu][%lu]=%f\n",i,j, dom.Rho[0][1][i][j]);

            }

            if (!dom.IsSolid[1][dom.Ndim(0) - 2][i][j])
            {
                double rho = dom.Rho[1][dom.Ndim(0) - 2][i][j]; // 同上
                rhoo += rho;
                nfo++;
            }
        }
    }

  
    
    // if (nfb > 0)
    rhow /= nfb;
    //   if (nfo > 0)
    rhoo /= nfo;

    double Pc;
    double rho;
     size_t newSum = (nw+no);
    double newSr = (1.00*nw)/newSum;
    double a = M_PI / dat.ome;
    rho = dat.Head * ((1.0 / a) * (dat.time - a * (floor(dat.time / a) + 0.5)) * pow(-1.0, floor(dat.time / a)) + 0.5) + dat.Orig;
    Pc = (2.0 * (rho - dat.rho) + dom.Gmix * (rho * rho * 0.999 * 0.001 - (2.0 * dat.rho - rho) * (2.0 * dat.rho - rho) * 0.999 * 0.001)) / 3.0;
   // dat.oss_ss << dom.Time  << Util::_8s << "nw"<< Util::_8s << nw << Util::_8s << "no"<< Util::_8s << no << Util::_8s << "Pc" << Pc << Util::_8s << Util::_8s << "Sr = nw/(nw+no)"<< newSr<< Util::_8s << "newSum = (nw+no)"<< Util::_8s << newSum << Util::_8s << "Sr-newSr"<< Util::_8s << Sr-newSr << Util::_8s << "Ncells- Sf"<< dom.Ncells-Sf<< std::endl;
    
    dat.oss_ss << dom.Time << Util::_8s << rho << Util::_8s << rhoo << Util::_8s << rhow << Util::_8s << water << Util::_8s << oil << Util::_8s << Pc << Util::_8s << Sr << std::endl;
}

int main(int argc, char **argv)
try
{
    String filekey(argv[1]);
    String filename(filekey + ".inp");
    if (!Util::FileExists(filename))
        throw new Fatal("File <%s> not found", filename.CStr());
    std::ifstream infile(filename.CStr());
    size_t Nproc = 1;
    if (argc == 3)
        Nproc = atoi(argv[2]);

    String fileDEM;
    String fileLBM;
    bool Render = true;
    size_t N = 200;
    real Gs0 = -0.53;
    real Gs1 = -0.53;
    real Gmix = 2.0;
    double nu = 0.05;
    double dt = 1.0;
    double Tf = 10000.0;
    real dtOut = 50.0;
    real HeadStep = 1000.0;
    real rho = 200.0;
    real ome = 2.0;
    real Head = 500.0;
    real Orig = 54.0;
    size_t oct = 1;
    real DPx = 1.0;
    real DPy = 1.0;
    real DPz = 1.0;
    int outlimit = 1;
    size_t buffer = 1;
    {
        infile >> fileDEM;
        infile.ignore(200, '\n');
        infile >> fileLBM;
        infile.ignore(200, '\n');
        infile >> Render;
        infile.ignore(200, '\n');
        infile >> N;
        infile.ignore(200, '\n');
        infile >> Gs0;
        infile.ignore(200, '\n');
        infile >> Gs1;
        infile.ignore(200, '\n');
        infile >> Gmix;
        infile.ignore(200, '\n');
        infile >> nu;
        infile.ignore(200, '\n');
        infile >> dt;
        infile.ignore(200, '\n');
        infile >> Tf;
        infile.ignore(200, '\n');
        infile >> dtOut;
        infile.ignore(200, '\n');
        infile >> HeadStep;
        infile.ignore(200, '\n');
        infile >> rho;
        infile.ignore(200, '\n');
        infile >> ome;
        infile.ignore(200, '\n');
        infile >> Head;
        infile.ignore(200, '\n');
        infile >> Orig;
        infile.ignore(200, '\n');
        infile >> oct;
        infile.ignore(200, '\n');
        infile >> DPx;
        infile.ignore(200, '\n');
        infile >> DPy;
        infile.ignore(200, '\n');
        infile >> DPz;
        infile.ignore(200, '\n');
        infile >> outlimit;
        infile.ignore(200, '\n');
        infile >> buffer;
        infile.ignore(200, '\n');
    }
    Array<real> nua(2);
    nua[0] = nu;
    nua[1] = nu;

    DEM::Domain DemDom;
    DemDom.Load(fileDEM.CStr());
    Array<int> idx(6);
    idx = -2, -3, -4, -5, -6, -7;
    DemDom.DelParticles(idx);
    Vec3_t Xmin, Xmax;
    DemDom.BoundingBox(Xmin, Xmax);
    int bound = outlimit;
    real dx = (Xmax(0) - Xmin(0)) / (N - 2 * bound);
    size_t Ny = (Xmax(1) - Xmin(1)) / dx + 2 * bound;
    size_t Nz = (Xmax(2) - Xmin(2)) / dx + 2 * bound;
    DemDom.Center(0.5 * (Xmax - Xmin) + Vec3_t(bound * dx, bound * dx, bound * dx));

    FLBM::Domain Dom(D3Q15, nua, iVec3_t(N, Ny, Nz), 1.0, 1.0);
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < Ny; y++)
        {
            for (int z = 0; z < Nz; z++)
            {
                Vec3_t pos((x + 0.5) * dx, (y + 0.5) * dx, (z + 0.5) * dx);
                for (DEM::Particle *P : DemDom.Particles)
                {
                    if (P->IsInsideAlt(pos))
                    {
                        Dom.IsSolid[0][x][y][z] = true;
                        Dom.IsSolid[1][x][y][z] = true;
                    }
                }
            }
        }
    }
    UserData dat;
    Dom.UserData = &dat;

    dat.Tf = Tf;
    dat.ome = 2 * M_PI * ome / Tf;
    dat.Orig = Orig;
    dat.dtOut = HeadStep;
    dat.time = 0.0;
    dat.rho = rho;
    dat.Head = Head;
    dat.Dp = Vec3_t(DPx, DPy, DPz);
    dat.Dp /= norm(dat.Dp);
    dat.block = oct;

    Dom.G[0] = 0.0;
    Dom.Gs[0] = Gs0;
    Dom.G[1] = 0.0;
    Dom.Gs[1] = Gs1;
    Dom.Gmix = Gmix;
    dat.Xmin.resize(2);
    dat.Xmax.resize(2);
    dat.Ymin.resize(2);
    dat.Ymax.resize(2);
    dat.Zmin.resize(2);
    dat.Zmax.resize(2);
    // The 6 faces (x,y) (x,z) (y,z) of the cube correspond to 2 of each group


    for (int i = 0; i < N; i++)
    {
        Dom.IsSolid[0][i][0][0] = true;
        Dom.IsSolid[0][i][Ny - 1][0] = true;
        Dom.IsSolid[0][i][0][Nz - 1] = true;
        Dom.IsSolid[0][i][Ny - 1][Nz - 1] = true;
        Dom.IsSolid[1][i][0][0] = true;
        Dom.IsSolid[1][i][Ny - 1][0] = true;
        Dom.IsSolid[1][i][0][Nz - 1] = true;
        Dom.IsSolid[1][i][Ny - 1][Nz - 1] = true;
    }

    for (int i = 0; i < Ny; i++)
    {
        Dom.IsSolid[0][0][i][0] = true;
        Dom.IsSolid[0][N - 1][i][0] = true;
        Dom.IsSolid[0][0][i][Nz - 1] = true;
        Dom.IsSolid[0][N - 1][i][Nz - 1] = true;
        Dom.IsSolid[1][0][i][0] = true;
        Dom.IsSolid[1][N - 1][i][0] = true;
        Dom.IsSolid[1][0][i][Nz - 1] = true;
        Dom.IsSolid[1][N - 1][i][Nz - 1] = true;
    }

    for (int i = 0; i < Nz; i++)
    {
        Dom.IsSolid[0][0][0][i] = true;
        Dom.IsSolid[0][N - 1][0][i] = true;
        Dom.IsSolid[0][0][Ny - 1][i] = true;
        Dom.IsSolid[0][N - 1][Ny - 1][i] = true;
        Dom.IsSolid[1][0][0][i] = true;
        Dom.IsSolid[1][N - 1][0][i] = true;
        Dom.IsSolid[1][0][Ny - 1][i] = true;
        Dom.IsSolid[1][N - 1][Ny - 1][i] = true;
    }

    bound = buffer;
    for (size_t ix = 0; ix < Dom.Ndim(0); ix++)
    {
        for (size_t iy = 0; iy < Dom.Ndim(1); iy++)
        {
            for (size_t iz = 0; iz < Dom.Ndim(2); iz++)
            {
                size_t il0 = 0;
                size_t il1 = 1;
                iVec3_t idx(ix, iy, iz);

                Dom.Initialize(il0, idx, 0.001 * rho, OrthoSys::O);
                Dom.Initialize(il1, idx, 0.999 * rho, OrthoSys::O);
                // if (oct < 2) // inp=>  oct=1

                {
                    //if ((dat.Dp(0) > 1.0e-12) && (ix < Dom.Ndim(0) / bound)) // 100     //change to 2
                    if ((dat.Dp(0) > 1.0e-12) && (ix < bound))
                    {
                        Dom.Initialize(il0, idx, 0.999 * rho, OrthoSys::O);
                        Dom.Initialize(il1, idx, 0.001 * rho, OrthoSys::O);
                    }
                   // if ((dat.Dp(1) > 1.0e-12) && (iy < Dom.Ndim(1) / bound)) // 90
                    if ((dat.Dp(1) > 1.0e-12) && (iy < bound)) // 90
                    {
                        Dom.Initialize(il0, idx, 0.999 * rho, OrthoSys::O);
                        Dom.Initialize(il1, idx, 0.001 * rho, OrthoSys::O);
                    }
                    //if ((dat.Dp(2) > 1.0e-12) && (iz < Dom.Ndim(2) / bound)) // 93
                     if ((dat.Dp(2) > 1.0e-12) && (iz < bound)) // 93
                    {
                        Dom.Initialize(il0, idx, 0.999 * rho, OrthoSys::O);
                        Dom.Initialize(il1, idx, 0.001 * rho, OrthoSys::O);
                    }
                    /*

                                        std::string output;
                                        for (int k = 0; k < 15; k++)
                                        {
                                            output += "F[1][ix][iy][iz][" + std::to_string(k) + "] = " + std::to_string(Dom.F[1][0][0][0][k]) + "\n";
                                        }

                                        printf("%s", output.c_str());
                    */
                }
            }
        }
    }

    // Dom.WriteXDMF("wreten");

    String fs;
    fs.Printf("water_retention.res");
    dat.oss_ss.open(fs.CStr(), std::ios::out);
    dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "PDen" << Util::_8s << "Rhow" << Util::_8s << "Rhoo" << Util::_8s << "Water" << Util::_8s << "Oil" << Util::_8s << "Pc" << Util::_8s << "Sr" << std::endl;
    Dom.Solve(Tf, dtOut, Setup, Report, filekey.CStr(), Render, Nproc);
    //Dom.Solve(Tf, dtOut, NULL, Report, filekey.CStr(), Render, Nproc);
    dat.oss_ss.close();
}
MECHSYS_CATCH
