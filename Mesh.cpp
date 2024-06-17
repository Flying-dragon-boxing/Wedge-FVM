#include <Eigen/Dense>
#include <cmath>
#include <fstream>

#include <cnpy.h>
#include <iostream>
#include <ostream>

#include "Mesh.h"

constexpr double GAMMA = 1.4;

Driver::Driver(Eigen::MatrixXd &x, Eigen::MatrixXd &y)
{
    this->x = x;
    this->y = y;
    this->nx = x.cols()-1;
    this->ny = x.rows()-1;

    // allocate memory for rho, u, v, p, note: two ghost cells are added to each side
    this->rho = Eigen::MatrixXd::Zero(ny+4, nx+4);
    this->u = Eigen::MatrixXd::Zero(ny+4, nx+4);
    this->v = Eigen::MatrixXd::Zero(ny+4, nx+4);
    this->p = Eigen::MatrixXd::Zero(ny+4, nx+4);

    this->fluxes.resize(nx+4);
    for (int i = 0; i < nx+4; i++)
    {
        this->fluxes[i].resize(ny+4);
    }

}

Driver::Driver()
{
    std::ifstream param("Mesh/param.ini");
    int nx, ny;
    param >> nx >> ny;

    // read mesh from file
    cnpy::NpyArray arr_x = cnpy::npy_load("Mesh/x.npy");
    cnpy::NpyArray arr_y = cnpy::npy_load("Mesh/y.npy");

    // convert to Eigen::MatrixXd
    x = Eigen::Map<Eigen::MatrixXd>(arr_x.data<double>(), ny, nx).transpose();
    y = Eigen::Map<Eigen::MatrixXd>(arr_y.data<double>(), ny, nx).transpose();

    this->ny = x.cols()-1;
    this->nx = x.rows()-1;

    nx -= 1;
    ny -= 1;

    // allocate memory for rho, u, v, p, note: two ghost cells are added to each side
    this->rho = Eigen::MatrixXd::Zero(nx+4, ny+4);
    this->u = Eigen::MatrixXd::Zero(nx+4, ny+4);
    this->v = Eigen::MatrixXd::Zero(nx+4, ny+4);
    this->p = Eigen::MatrixXd::Zero(nx+4, ny+4);

    this->fluxes.resize(nx+4);
    for (int i = 0; i < nx+4; i++)
    {
        this->fluxes[i].resize(ny+4);
    }
}

void Driver::apply_initial_conditions()
{
    #pragma omp parallel for schedule(static)
    for (int i = 2; i < nx + 2; i++)
    {
        for (int j = 2; j < ny + 2; j++)
        {
            // rho(i, j) = 2.6662283417282984;
            // u(i, j) = 0.7501232991558267;
            // v(i, j) = 0.0;
            // p(i, j) = 3.2142504859554784;

            // rho(i, j) = 1.1846616239476357;
            // u(i, j) = 686.47;
            // v(i, j) = 0.0;
            // p(i, j) = 99719;

            rho(i, j) = 1.0;
            u(i, j) = 2.0;
            v(i, j) = 0.0;
            p(i, j) = 0.714497;
            
        }
    }
}

void Driver::apply_boundary_conditions()
{
    Driver &mesh = *this;

    // left: inlet, right: outlet, top: outlet, bottom: wall
    // update only ghost cells

    // left
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < mesh.ny + 4; j++) 
    {
        mesh.rho(0, j) = 1.0;
        mesh.u(0, j) = 2.0;
        mesh.v(0, j) = 0.0;
        mesh.p(0, j) = 0.714497;

        mesh.rho(1, j) = 1.0;
        mesh.u(1, j) = 2.0;
        mesh.v(1, j) = 0.0;
        mesh.p(1, j) = 0.714497;

        // mesh.rho(0, j) = 2.6662283417282984;
        // mesh.u(0, j) = 0.7501232991558267;
        // mesh.v(0, j) = 0.0;
        // mesh.p(0, j) = 3.2142504859554784;

        // mesh.rho(1, j) = 2.6662283417282984;
        // mesh.u(1, j) = 0.7501232991558267;
        // mesh.v(1, j) = 0.0;
        // mesh.p(1, j) = 3.2142504859554784;
    }

    // right
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < mesh.ny + 4; j++)
    {
        mesh.rho(mesh.nx + 2, j) = mesh.rho(mesh.nx + 1, j);
        mesh.u(mesh.nx + 2, j) = mesh.u(mesh.nx + 1, j);
        mesh.v(mesh.nx + 2, j) = mesh.v(mesh.nx + 1, j);
        mesh.p(mesh.nx + 2, j) = mesh.p(mesh.nx + 1, j);

        mesh.rho(mesh.nx + 3, j) = mesh.rho(mesh.nx + 1, j);
        mesh.u(mesh.nx + 3, j) = mesh.u(mesh.nx + 1, j);
        mesh.v(mesh.nx + 3, j) = mesh.v(mesh.nx + 1, j);
        mesh.p(mesh.nx + 3, j) = mesh.p(mesh.nx + 1, j);
    }

    // top
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < mesh.nx + 4; i++)
    {
        mesh.rho(i, mesh.ny + 2) = mesh.rho(i, mesh.ny + 1);
        mesh.u(i, mesh.ny + 2) = mesh.u(i, mesh.ny + 1);
        mesh.v(i, mesh.ny + 2) = mesh.v(i, mesh.ny + 1);
        mesh.p(i, mesh.ny + 2) = mesh.p(i, mesh.ny + 1);

        mesh.rho(i, mesh.ny + 3) = mesh.rho(i, mesh.ny + 1);
        mesh.u(i, mesh.ny + 3) = mesh.u(i, mesh.ny + 1);
        mesh.v(i, mesh.ny + 3) = mesh.v(i, mesh.ny + 1);
        mesh.p(i, mesh.ny + 3) = mesh.p(i, mesh.ny + 1);
    }

    // bottom
    // #pragma omp parallel for schedule(static)
    // for (int i = 0; i < 22; i++)
    // {
    //     mesh.rho(i, 0) = mesh.rho(i, 2);
    //     mesh.u(i, 0) = mesh.u(i, 2);
    //     mesh.v(i, 0) = mesh.v(i, 2);
    //     mesh.p(i, 0) = mesh.p(i, 2);

    //     mesh.rho(i, 1) = mesh.rho(i, 2);
    //     mesh.u(i, 1) = mesh.u(i, 2);
    //     mesh.v(i, 1) = mesh.v(i, 2);
    //     mesh.p(i, 1) = mesh.p(i, 2);
    // }

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < mesh.nx + 4; i++)
    {
        mesh.rho(i, 0) = mesh.rho(i, 2);
        if ((i > 2 && i < mesh.nx + 2 && mesh.x(i-1, 0) > 1))
        {        
            // mesh.u(i, 0) = -mesh.u(i, 2);
            // mesh.v(i, 0) = -mesh.v(i, 2);
            Vec2 v0 = {mesh.u(i, 2), mesh.v(i, 2)};
            Vec2 tangent = {x(i-1, 0) - x(i-2, 0), y(i-1, 0) - y(i-2, 0)};
            tangent = tangent.normalize();
            Vec2 normal = {tangent.y, -tangent.x};
            double v0_n = v0.dot(normal);
            double v0_t = v0.dot(tangent);

            double v_t = v0_t;
            double v_n = -v0_n;
            Vec2 v_new = {v_n * normal.x + v_t * tangent.x, v_n * normal.y + v_t * tangent.y};
            mesh.u(i, 0) = v_new.x;
            mesh.v(i, 0) = v_new.y;
        }
        else 
        {
            mesh.u(i, 0) = 2;//-mesh.u(i, 2) + 4;
            mesh.v(i, 0) = 0;//-mesh.v(i, 2);
        }
        mesh.p(i, 0) = mesh.p(i, 2);

        mesh.rho(i, 1) = mesh.rho(i, 2);
        if ((i > 2 && i < mesh.nx + 2 && mesh.x(i-1, 0) > 1))
        {        
            // mesh.u(i, 1) = -mesh.u(i, 2);
            // mesh.v(i, 1) = -mesh.v(i, 2);
            Vec2 v0 = {mesh.u(i, 2), mesh.v(i, 2)};
            Vec2 tangent = {x(i-1, 0) - x(i-2, 0), y(i-1, 0) - y(i-2, 0)};
            tangent = tangent.normalize();
            Vec2 normal = {tangent.y, -tangent.x};
            double v0_n = v0.dot(normal);
            double v0_t = v0.dot(tangent);

            double v_t = v0_t;
            double v_n = -v0_n;
            Vec2 v_new = {v_n * normal.x + v_t * tangent.x, v_n * normal.y + v_t * tangent.y};
            mesh.u(i, 1) = v_new.x;
            mesh.v(i, 1) = v_new.y;
        }
        else 
        {
            mesh.u(i, 1) = 2;//-mesh.u(i, 2) + 4;
            mesh.v(i, 1) = 0;//-mesh.v(i, 2);
        }
        mesh.p(i, 1) = mesh.p(i, 2);
    }
}

bool Driver::solve(double dt)
{
    Eigen::MatrixXd rho_inc = Eigen::MatrixXd::Zero(nx+4, ny+4);
    Eigen::MatrixXd u_inc = Eigen::MatrixXd::Zero(nx+4, ny+4);
    Eigen::MatrixXd v_inc = Eigen::MatrixXd::Zero(nx+4, ny+4);
    Eigen::MatrixXd p_inc = Eigen::MatrixXd::Zero(nx+4, ny+4);

    bool negative_pressure = false;

    // for each inner cell
    #pragma omp parallel for schedule(static) collapse(2) reduction(||:negative_pressure)
    for (int i = 2; i < nx + 2; i++)
    {
        for (int j = 2; j < ny + 2; j++)
        {
            // unit points
            Vec2 nw = {x(i-2, j-1), y(i-2, j-1)};
            Vec2 ne = {x(i-1, j-1), y(i-1, j-1)};
            Vec2 sw = {x(i-2, j-2), y(i-2, j-2)};
            Vec2 se = {x(i-1, j-2), y(i-1, j-2)};

            // unit volume
            double vol = ((ne - nw).dot(se - nw) + (se - nw).dot(sw - nw)) / 2.0;

            // // north
            // {
            //     double rho_n = 0.5 * (rho(i, j) + rho(i, j+1));
            //     double u_n = 0.5 * (u(i, j) + u(i, j+1));
            //     double v_n = 0.5 * (v(i, j) + v(i, j+1));
            //     double p_n = 0.5 * (p(i, j) + p(i, j+1));

            //     // Vec2 n_bnd = {x(i-1, j-1) - x(i-2, j-1), y(i-1, j-1) - y(i-2, j-1)};
            //     Vec2 n_bnd = ne - nw;
            //     double len = n_bnd.norm();
            //     n_bnd = n_bnd.normalize();
            //     Vec2 n_normal = {n_bnd.y, -n_bnd.x};

            //     Vec2 u_bnd = {u_n, v_n};
            //     double u_normal = u_bnd.dot(n_normal);
            //     double u_tangent = u_bnd.dot(n_bnd);

            //     // // reconstruction, eigenvalues
            //     // double a = std::sqrt(1.4 * p_n / rho_n);
            //     // double u_n_reconstruct = 1;

            //     // riemann solver

            //     // roe average
            //     double rho_roe = std::sqrt(rho(i, j) * rho(i, j+1));
            //     double u_roe = (std::sqrt(rho(i, j)) * u(i, j) + std::sqrt(rho(i, j+1)) * u(i, j+1)) / (std::sqrt(rho(i, j)) + std::sqrt(rho(i, j+1)));
                
            // }

            double rho_u_inc = 0.0;
            double rho_v_inc = 0.0;
            double e_inc = 0.0;

            // north
            {
                Vec2 n_bnd = ne - nw;
                double len = n_bnd.norm();
                n_bnd = n_bnd.normalize();
                Vec2 n_normal = {-n_bnd.y, n_bnd.x};

                Vec2 u_bnd = {u(i, j), v(i, j)};
                double u_normal = u_bnd.dot(n_normal);
                double u_tangent = u_bnd.dot(n_bnd);

                rho_inc(i, j) -= dt / vol * fluxes[i][j].n[0] * len;
                double rho_u_inc_tangent = - dt / vol * fluxes[i][j].n[1] * len;
                double rho_u_inc_normal = - dt / vol * fluxes[i][j].n[2] * len;
                e_inc -= dt / vol * fluxes[i][j].n[3] * len;

                // transform the velocity increment to the cartesian coordinate
                rho_u_inc += rho_u_inc_normal * n_normal.x + rho_u_inc_tangent * n_bnd.x;
                rho_v_inc += rho_u_inc_normal * n_normal.y + rho_u_inc_tangent * n_bnd.y;
            }

            // east
            {
                Vec2 e_bnd = se - ne;
                double len = e_bnd.norm();
                e_bnd = e_bnd.normalize();
                Vec2 e_normal = {-e_bnd.y, e_bnd.x};

                Vec2 u_bnd = {u(i, j), v(i, j)};
                double u_normal = u_bnd.dot(e_normal);
                double u_tangent = u_bnd.dot(e_bnd);

                rho_inc(i, j) -= dt / vol * fluxes[i][j].e[0] * len;
                double rho_u_inc_tangent = - dt / vol * fluxes[i][j].e[1] * len;
                double rho_u_inc_normal = - dt / vol * fluxes[i][j].e[2] * len;
                e_inc -= dt / vol * fluxes[i][j].e[3] * len;

                // transform the velocity increment to the cartesian coordinate
                rho_u_inc += rho_u_inc_normal * e_normal.x + rho_u_inc_tangent * e_bnd.x;
                rho_v_inc += rho_u_inc_normal * e_normal.y + rho_u_inc_tangent * e_bnd.y;
            }

            // south
            {
                Vec2 s_bnd = sw - se;
                double len = s_bnd.norm();
                s_bnd = s_bnd.normalize();
                Vec2 s_normal = {-s_bnd.y, s_bnd.x};

                Vec2 u_bnd = {u(i, j), v(i, j)};
                double u_normal = u_bnd.dot(s_normal);
                double u_tangent = u_bnd.dot(s_bnd);

                rho_inc(i, j) -= dt / vol * fluxes[i][j].s[0] * len;
                double rho_u_inc_tangent = - dt / vol * fluxes[i][j].s[1] * len;
                double rho_u_inc_normal = - dt / vol * fluxes[i][j].s[2] * len;
                e_inc -= dt / vol * fluxes[i][j].s[3] * len;

                // transform the velocity increment to the cartesian coordinate
                rho_u_inc += rho_u_inc_normal * s_normal.x + rho_u_inc_tangent * s_bnd.x;
                rho_v_inc += rho_u_inc_normal * s_normal.y + rho_u_inc_tangent * s_bnd.y;
            }

            // west
            {
                Vec2 w_bnd = nw - sw;
                double len = w_bnd.norm();
                w_bnd = w_bnd.normalize();
                Vec2 w_normal = {-w_bnd.y, w_bnd.x};

                Vec2 u_bnd = {u(i, j), v(i, j)};
                double u_normal = u_bnd.dot(w_normal);
                double u_tangent = u_bnd.dot(w_bnd);

                rho_inc(i, j) -= dt / vol * fluxes[i][j].w[0] * len;
                double rho_u_inc_tangent = - dt / vol * fluxes[i][j].w[1] * len;
                double rho_u_inc_normal = - dt / vol * fluxes[i][j].w[2] * len;
                e_inc -= dt / vol * fluxes[i][j].w[3] * len;

                // transform the velocity increment to the cartesian coordinate
                rho_u_inc += rho_u_inc_normal * w_normal.x + rho_u_inc_tangent * w_bnd.x;
                rho_v_inc += rho_u_inc_normal * w_normal.y + rho_u_inc_tangent * w_bnd.y;
            }

            if (std::isnan(rho_u_inc))
            {
                rho_u_inc = 0;
                // std::cout << "nan" << std::endl;
                // #pragma omp critical
                //     return true;
            }
            if (std::isnan(rho_v_inc))
            {
                rho_v_inc = 0;
                // std::cout << "nan" << std::endl;
                // #pragma omp critical
                //     return true;
            }
            if (std::isnan(e_inc))
            {
                e_inc = 0;
                // std::cout << "nan" << std::endl;
                // #pragma omp critical
                //     return true;
            }
            if (std::isnan(rho_inc(i, j)))
            {
                rho_inc(i, j) = 0;
                // std::cout << "nan" << std::endl;
                // #pragma omp critical
                //     return true;
            }

            u_inc(i, j) = rho_u_inc / (rho(i, j) + rho_inc(i, j));
            v_inc(i, j) = rho_v_inc / (rho(i, j) + rho_inc(i, j));
            double u_new = u(i, j) + u_inc(i, j);
            double v_new = v(i, j) + v_inc(i, j);
            double rho_new = rho(i, j) + rho_inc(i, j);
            double e = p(i, j) / (GAMMA - 1) + 0.5 * rho(i, j) * (u(i, j) * u(i, j) + v(i, j) * v(i, j));
            double e_new = e + e_inc;
            // double p_new = (GAMMA - 1) * (e_new - 0.5 * rho(i, j) * (u(i, j) * u(i, j) + v(i, j) * v(i, j)));
            double p_new = (GAMMA - 1) * (e_new - 0.5 * rho(i,j) * (u_new * u_new + v_new * v_new));
            if (p_new < 0)
            {
                std::cout << "negative pressure" << std::endl;
                negative_pressure = true;
            }
            p_inc(i, j) = p_new - p(i, j);
            // // std::cout << "p_inc(" << i<< ", " << j << "): " << p_inc(i, j) << std::endl;
            p(i, j) = p_new;
            // if (i == 60 && j == 2)
            // {
            //     std::cout << "e_inc: " << e_inc << std::endl;
            //     std::cout << "p_inc: " << p_inc(i, j) << std::endl;
            // }
        }

    }

    rho += rho_inc;
    u += u_inc;
    v += v_inc;
    bool converge = (rho_inc.cwiseAbs().maxCoeff() < threshold && u_inc.cwiseAbs().maxCoeff() < threshold && v_inc.cwiseAbs().maxCoeff() < threshold ) || negative_pressure;
    return converge;

    
}
void Driver::add_viscosity()
{
    double ratio = 0.2;
    // we do this by mixing every cell's value with the average of its neighbors
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 2; i < nx+2; i++)
    {
        for (int j = 2; j < 10; j++)
        {
            rho(i, j) = 0.25 * ratio * (rho(i-1, j) + rho(i+1, j) + rho(i, j-1) + rho(i, j+1)) + (1 - ratio) * rho(i, j);
            u(i, j) = 0.25 * ratio * (u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1)) + (1 - ratio) * u(i, j);
            v(i, j) = 0.25 * ratio * (v(i-1, j) + v(i+1, j) + v(i, j-1) + v(i, j+1)) + (1 - ratio) * v(i, j);
            p(i, j) = 0.25 * ratio * (p(i-1, j) + p(i+1, j) + p(i, j-1) + p(i, j+1)) + (1 - ratio) * p(i, j);
        }
    }
}

void Driver::save()
{
    cnpy::npy_save("rho.npy", rho.data(), {static_cast<unsigned long>(ny+4), static_cast<unsigned long>(nx+4)}, "w");
    cnpy::npy_save("u.npy", u.data(), {static_cast<unsigned long>(ny+4), static_cast<unsigned long>(nx+4)}, "w");
    cnpy::npy_save("v.npy", v.data(), {static_cast<unsigned long>(ny+4), static_cast<unsigned long>(nx+4)}, "w");
    cnpy::npy_save("p.npy", p.data(), {static_cast<unsigned long>(ny+4), static_cast<unsigned long>(nx+4)}, "w");
}

double Vec2::dot(Vec2 v)
{
    return this->x * v.x + this->y * v.y;
}

Vec2 Vec2::operator+(Vec2 v)
{
    return {this->x + v.x, this->y + v.y};
}

Vec2 Vec2::operator-(Vec2 v)
{
    return {this->x - v.x, this->y - v.y};
}

Vec2 Vec2::operator*(double s)
{
    return {this->x * s, this->y * s};
}

Vec2 Vec2::operator/(double s)
{
    return {this->x / s, this->y / s};
}

double Vec2::norm()
{
    return std::sqrt(this->x * this->x + this->y * this->y);
}

Vec2 Vec2::normalize()
{
    assert (this->norm() > 1e-10);
    return *this / this->norm();
}
