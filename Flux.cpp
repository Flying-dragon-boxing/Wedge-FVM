#include <Eigen/Dense>
#include <cmath>
#include <fstream>

#include <cnpy.h>
#include <iostream>
#include <ostream>

#include "Mesh.h"

constexpr double GAMMA = 1.4;

inline double minmod(const double &a, const double &b)
{
    bool a_b_same_sign = a * b > 0;
    if (a_b_same_sign)
    {
        return std::abs(a) < std::abs(b) ? a : b;
    }
    else
    {
        return 0.0;
    }

}

inline double minmod(const double &a, const double &b, const double &c)
{
    return minmod(a, minmod(b, c));
}

inline double Driver::entropy_abs(const double &a, const double &dx, const double &dt)
{
    double new_a = a * dt / dx;
    if (new_a > epsilon)
    {
        return new_a * dx / dt;
    }
    else if (new_a < -epsilon)
    {
        return -new_a * dx / dt;
    }
    else
    {
        // return (new_a * new_a + epsilon * epsilon) / 2 / epsilon * dx / dt;
        return epsilon * dx / dt;
    }
}

Flux Driver::reconstruct_v_l(const State &mm, const State &m, const State &p, const State &pp, Vec2 normal)
{
    // m: minus, p: plus
    // Here shows how the units are arranged
    // mm | m | p | pp
    // normal is ->
    double u_mm = mm[1], u_m = m[1], u_p = p[1], u_pp = pp[1];
    double v_mm = mm[2], v_m = m[2], v_p = p[2], v_pp = pp[2];
    double a_mm = std::sqrt(GAMMA * mm[3] / mm[0]),
            a_m = std::sqrt(GAMMA * m[3] / m[0]),
            a_p = std::sqrt(GAMMA * p[3] / p[0]),
           a_pp = std::sqrt(GAMMA * pp[3] / pp[0]);

    // convert to normal coordinate
    double u_mm_n = u_mm * normal.x + v_mm * normal.y, u_mm_t = u_mm * normal.y - v_mm * normal.x;
    double u_m_n = u_m * normal.x + v_m * normal.y, u_m_t = u_m * normal.y - v_m * normal.x;
    double u_p_n = u_p * normal.x + v_p * normal.y, u_p_t = u_p * normal.y - v_p * normal.x;
    double u_pp_n = u_pp * normal.x + v_pp * normal.y, u_pp_t = u_pp * normal.y - v_pp * normal.x;

    // apply original variable reconstruction
    double rho_m = m[0];// + 0.5 * minmod(m[0] - mm[0], p[0] - m[0]);
    double rho_p = p[0];// - 0.5 * minmod( p[0] - m[0], pp[0] - p[0]);
    u_m_n = u_m_n + 0.5 * minmod(u_m_n - u_mm_n, u_p_n - u_m_n);
    u_p_n = u_p_n - 0.5 * minmod(u_p_n - u_m_n, u_pp_n - u_p_n);
    u_m_t = u_m_t + 0.5 * minmod(u_m_t - u_mm_t, u_p_t - u_m_t);
    u_p_t = u_p_t - 0.5 * minmod(u_p_t - u_m_t, u_pp_t - u_p_t);
    double pressure_m = m[3];// + 0.5 * minmod(m[3] - mm[3], p[3] - m[3]);
    double pressure_p = p[3];// - 0.5 * minmod(p[3] - m[3], pp[3] - p[3]);

    a_m = std::sqrt(GAMMA * pressure_m / rho_m);
    a_p = std::sqrt(GAMMA * pressure_p / rho_p);

    // experimental
    double a_mid = (rho_m * a_m + rho_p * a_p) / (rho_m + rho_p);
    a_m = a_mid; a_p = a_mid;

    double mach_normal_m = u_m_n / a_m, mach_normal_p = u_p_n / a_p;

    auto soft_pos = [&](const double &a) {
      if (a > 1) {
        return a;
      } else if (a < -1) {
        return 0.0;
      } else {
        return 0.25 * (a + 1) * (a + 1);
      }
    };

    auto soft_neg = [&](const double &a) {
      if (a < -1) {
        return a;
      } else if (a > 1) {
        return 0.0;
      } else {
        return -0.25 * (a - 1) * (a - 1);
      }
    };

    mach_normal_m = soft_pos(mach_normal_m);
    mach_normal_p = soft_neg(mach_normal_p);
    double mach_normal = mach_normal_m + mach_normal_p;

    if (mach_normal < -1)
    {
        return {rho_p * u_p_n, rho_p * u_p_n * u_p_t, rho_p * u_p_n * u_p_n + pressure_p, u_p_n * (pressure_p + 1/(GAMMA - 1) * pressure_p + 0.5 * rho_p * (u_p_n * u_p_n + u_p_t * u_p_t))};
    }
    else if (mach_normal > 1)
    {
        return {rho_m * u_m_n, rho_m * u_m_n * u_m_t, rho_m * u_m_n * u_m_n + pressure_m, u_m_n * (pressure_m + 1/(GAMMA - 1) * pressure_m + 0.5 * rho_m * (u_m_n * u_m_n + u_m_t * u_m_t))};
    }

    auto square = [](const double &a) { return a * a; };

    double f_m = rho_m * a_m * mach_normal_m;
    double f_p = rho_p * a_p * mach_normal_p;

    
    // u_m_n = u_n, u_p_n = u_n;
    // double u_n = mach_normal * a_mid;


    Flux flux_m = {
        f_m,
        f_m * u_m_t,
        f_m * (u_m_n + (-u_m_n + 2 * a_m) / GAMMA),
        f_m * (square((GAMMA-1)*u_m_n + 2 * a_m) / 2 / (GAMMA * GAMMA - 1) + (u_m_n * u_m_n + u_m_t * u_m_t - u_m_n * u_m_n) / 2)
    };

    Flux flux_p = {
        f_p,
        f_p * u_p_t,
        f_p * (u_m_n + (-u_m_n - 2 * a_p) / GAMMA),
        f_p * (square((GAMMA-1)*u_p_n - 2 * a_p) / 2 / (GAMMA * GAMMA - 1) + (u_p_n * u_p_n + u_p_t * u_p_t - u_p_n * u_p_n) / 2)
    };


    Flux flux = {
        flux_m[0] + flux_p[0],
        flux_m[1] + flux_p[1],
        flux_m[2] + flux_p[2],
        flux_m[3] + flux_p[3]
    };
    return flux;

}

Flux Driver::reconstruct_roe(const State &mm, const State &m, const State &p, const State &pp, Vec2 normal)
{
    // m: minus, p: plus
    // Here shows how the units are arranged
    // mm | m | p | pp
    // normal is ->
    double u_mm = mm[1], u_m = m[1], u_p = p[1], u_pp = pp[1];
    double v_mm = mm[2], v_m = m[2], v_p = p[2], v_pp = pp[2];
    double a_mm = std::sqrt(GAMMA * mm[3] / mm[0]),
            a_m = std::sqrt(GAMMA * m[3] / m[0]),
            a_p = std::sqrt(GAMMA * p[3] / p[0]),
           a_pp = std::sqrt(GAMMA * pp[3] / pp[0]);

    // convert to normal coordinate
    double u_mm_n = u_mm * normal.x + v_mm * normal.y, u_mm_t = u_mm * normal.y - v_mm * normal.x;
    double u_m_n = u_m * normal.x + v_m * normal.y, u_m_t = u_m * normal.y - v_m * normal.x;
    double u_p_n = u_p * normal.x + v_p * normal.y, u_p_t = u_p * normal.y - v_p * normal.x;
    double u_pp_n = u_pp * normal.x + v_pp * normal.y, u_pp_t = u_pp * normal.y - v_pp * normal.x;

    double rho_m = m[0];
    double rho_p = p[0];
    double pressure_m = m[3];
    double pressure_p = p[3];

    // these are reconstructed velocities for eigenvalues
    double u_m_n_reconstruct = u_m_n + 0.5 * minmod(u_m_n - u_mm_n, u_p_n - u_m_n);
    double u_p_n_reconstruct = u_p_n - 0.5 * minmod(u_p_n - u_m_n, u_pp_n - u_p_n);
    double u_m_t_reconstruct = u_m_t + 0.5 * minmod(u_m_t - u_mm_t, u_p_t - u_m_t);
    double u_p_t_reconstruct = u_p_t - 0.5 * minmod(u_p_t - u_m_t, u_pp_t - u_p_t);
    double a_m_reconstruct = a_m + 0.5 * minmod(a_m - a_mm, a_p - a_m);
    double a_p_reconstruct = a_p - 0.5 * minmod(a_p - a_m, a_pp - a_p);
    double u_n_reconstruct = (std::sqrt(rho_m) * u_m_n_reconstruct + std::sqrt(rho_p) * u_p_n_reconstruct) / (std::sqrt(rho_m) + std::sqrt(rho_p));
    double u_t_reconstruct = (std::sqrt(rho_m) * u_m_t_reconstruct + std::sqrt(rho_p) * u_p_t_reconstruct) / (std::sqrt(rho_m) + std::sqrt(rho_p));
    double a_reconstruct = (std::sqrt(rho_m) * a_m_reconstruct + std::sqrt(rho_p) * a_p_reconstruct) / (std::sqrt(rho_m) + std::sqrt(rho_p));
    
    u_m_n = u_m_n_reconstruct;
    u_p_n = u_p_n_reconstruct;
    u_m_t = u_m_t_reconstruct;
    u_p_t = u_p_t_reconstruct;
    a_m = a_m_reconstruct;
    a_p = a_p_reconstruct;

    // Riemann solver Roe
    double /*rho_m = m[0], rho_p = p[0], */ rho_mid = std::sqrt(rho_m * rho_p);
    double u_n_mid = (std::sqrt(rho_m) * u_m_n + std::sqrt(rho_p) * u_p_n) / (std::sqrt(rho_m) + std::sqrt(rho_p));
    double u_t_mid = (std::sqrt(rho_m) * u_m_t + std::sqrt(rho_p) * u_p_t) / (std::sqrt(rho_m) + std::sqrt(rho_p));
    double h_m = GAMMA / (GAMMA - 1) * m[3] / rho_m + 0.5 * (u_m * u_m + v_m * v_m);
    double h_p = GAMMA / (GAMMA - 1) * p[3] / rho_m + 0.5 * (u_p * u_p + v_p * v_p);
    double h_mid = (std::sqrt(rho_m) * h_m + std::sqrt(rho_p) * h_p) / (std::sqrt(rho_m) + std::sqrt(rho_p));

    // double pressure_m = m[3];
    double e_m = 1/(GAMMA - 1) * pressure_m + 0.5 * rho_m * (u_m_n * u_m_n + u_m_t * u_m_t);
    // Flux flux = {rho * u_n, rho * u_n * u_t, rho * u_n * u_n + pressure, u_n * (pressure + e)};
    Flux flux_m = {rho_m * u_m_n, rho_m * u_m_n * u_m_t, rho_m * u_m_n * u_m_n + pressure_m, u_m_n * (pressure_m + e_m)};

    // double pressure_p = p[3];
    double e_p = 1/(GAMMA - 1) * pressure_p + 0.5 * rho_p * (u_p_n * u_p_n + u_p_t * u_p_t);
    // Flux flux = {rho * u_n, rho * u_n * u_t, rho * u_n * u_n + pressure, u_n * (pressure + e)};
    Flux flux_p = {rho_p * u_p_n, rho_p * u_p_n * u_p_t, rho_p * u_p_n * u_p_n + pressure_p, u_p_n * (pressure_p + e_p)};

    // ===========================================================
    // A naive implementation: average the fluxes on the both side
    // ===========================================================

    Flux flux = {0.5*(flux_m[0] + flux_p[0]), 0.5*(flux_m[1] + flux_p[1]), 0.5*(flux_m[2] + flux_p[2]), 0.5*(flux_m[3] + flux_p[3])};

    // ===========================================================
    // Roe's average: need another term in the flux
    // ===========================================================
    typedef Eigen::Matrix<double, 4, 4> Matrix4d;
    typedef Eigen::Matrix<double, 4, 1> Vector4d;
    Matrix4d A_abs = Matrix4d::Zero(), Lambda_abs = Matrix4d::Zero();
    auto e_abs = [&](const double &a) { return entropy_abs(a, std::abs(x(1, 0) - x(0, 0)), dt); };
    Matrix4d L = Matrix4d::Zero(), R = Matrix4d::Zero();

    // Private vars here
    {
        double u = u_t_mid, v = u_n_mid;
        double h = h_mid;
        double a = std::sqrt((GAMMA - 1) * (h - 0.5 * (u * u + v * v)));
        Lambda_abs(0, 0) = e_abs(u_n_reconstruct), Lambda_abs(1, 1) = e_abs(u_n_reconstruct), Lambda_abs(2, 2) = e_abs(u_n_reconstruct - a_reconstruct), Lambda_abs(3, 3) = e_abs(u_n_reconstruct + a_reconstruct);

        R << 0, 1, 1, 1,
             1, 0, u, u,
             0, v, v-a, v+a,
             u, 0.5*(v*v-u*u), h-v*a, h+v*a;

        double b2 = (GAMMA - 1) / a / a, b1 = 0.5 * b2 * (u * u + v * v);
        L << -b1*u, 1+b2*u*u, b2*u*v, -b2*u,
             1-b1, b2*u, b2*v, -b2,
             0.5*(b1+v/a), -0.5*b2*u, -0.5*(b2*v+1/a), 0.5*b2,
             0.5*(b1-v/a), -0.5*b2*u, -0.5*(b2*v-1/a), 0.5*b2;
    }

    A_abs = R * Lambda_abs * L;

    // Compute the flux
    Vector4d U_m;
    U_m << m[0], m[0] * u_m_t, m[0] * u_m_n, e_m;
    Vector4d U_p;
    U_p << p[0], p[0] * u_p_t, p[0] * u_p_n, e_p;
    Vector4d Flux_Roe_add = -0.5 * (A_abs * (U_p - U_m));
    flux = {flux[0] + Flux_Roe_add[0], flux[1] + Flux_Roe_add[1], flux[2] + Flux_Roe_add[2], flux[3] + Flux_Roe_add[3]};

    return flux;

}

Flux Driver::reconstruct_roe_boundary(const State &mm, const State &m, const State &p, const State &pp, Vec2 normal)
{
    // m: minus, p: plus
    // Here shows how the units are arranged
    // mm | m | p | pp
    // normal is ->
    double u_mm = mm[1], u_m = m[1], u_p = p[1], u_pp = pp[1];
    double v_mm = mm[2], v_m = m[2], v_p = p[2], v_pp = pp[2];
    double a_mm = std::sqrt(GAMMA * mm[3] / mm[0]),
            a_m = std::sqrt(GAMMA * m[3] / m[0]),
            a_p = std::sqrt(GAMMA * p[3] / p[0]),
           a_pp = std::sqrt(GAMMA * pp[3] / pp[0]);

    // convert to normal coordinate
    double u_mm_n = u_mm * normal.x + v_mm * normal.y, u_mm_t = u_mm * normal.y - v_mm * normal.x;
    double u_m_n = u_m * normal.x + v_m * normal.y, u_m_t = u_m * normal.y - v_m * normal.x;
    double u_p_n = u_p * normal.x + v_p * normal.y, u_p_t = u_p * normal.y - v_p * normal.x;
    double u_pp_n = u_pp * normal.x + v_pp * normal.y, u_pp_t = u_pp * normal.y - v_pp * normal.x;

    double rho_m = m[0];
    double rho_p = p[0];
    double pressure_m = m[3];
    double pressure_p = p[3];

    // these are reconstructed velocities for eigenvalues
    double u_m_n_reconstruct = u_m_n + 0.5 * minmod(u_m_n - u_mm_n, u_p_n - u_m_n);
    double u_p_n_reconstruct = u_p_n - 0.5 * minmod(u_p_n - u_m_n, u_pp_n - u_p_n);
    double u_m_t_reconstruct = u_m_t + 0.5 * minmod(u_m_t - u_mm_t, u_p_t - u_m_t);
    double u_p_t_reconstruct = u_p_t - 0.5 * minmod(u_p_t - u_m_t, u_pp_t - u_p_t);
    double a_m_reconstruct = a_m + 0.5 * minmod(a_m - a_mm, a_p - a_m);
    double a_p_reconstruct = a_p - 0.5 * minmod(a_p - a_m, a_pp - a_p);
    double u_n_reconstruct = (std::sqrt(rho_m) * u_m_n_reconstruct + std::sqrt(rho_p) * u_p_n_reconstruct) / (std::sqrt(rho_m) + std::sqrt(rho_p));
    double u_t_reconstruct = (std::sqrt(rho_m) * u_m_t_reconstruct + std::sqrt(rho_p) * u_p_t_reconstruct) / (std::sqrt(rho_m) + std::sqrt(rho_p));
    double a_reconstruct = (std::sqrt(rho_m) * a_m_reconstruct + std::sqrt(rho_p) * a_p_reconstruct) / (std::sqrt(rho_m) + std::sqrt(rho_p));

    // Riemann solver Roe
    double /*rho_m = m[0], rho_p = p[0], */ rho_mid = std::sqrt(rho_m * rho_p);
    double u_n_mid = (std::sqrt(rho_m) * u_m_n + std::sqrt(rho_p) * u_p_n) / (std::sqrt(rho_m) + std::sqrt(rho_p));
    double u_t_mid = (std::sqrt(rho_m) * u_m_t + std::sqrt(rho_p) * u_p_t) / (std::sqrt(rho_m) + std::sqrt(rho_p));
    double h_m = GAMMA / (GAMMA - 1) * m[3] / rho_m + 0.5 * (u_m * u_m + v_m * v_m);
    double h_p = GAMMA / (GAMMA - 1) * p[3] / rho_m + 0.5 * (u_p * u_p + v_p * v_p);
    double h_mid = (std::sqrt(rho_m) * h_m + std::sqrt(rho_p) * h_p) / (std::sqrt(rho_m) + std::sqrt(rho_p));

    // double pressure_m = m[3];
    double e_m = 1/(GAMMA - 1) * pressure_m + 0.5 * rho_m * (u_m_n * u_m_n + u_m_t * u_m_t);
    // Flux flux = {rho * u_n, rho * u_n * u_t, rho * u_n * u_n + pressure, u_n * (pressure + e)};
    Flux flux_m = {rho_m * u_m_n, rho_m * u_m_n * u_m_t, rho_m * u_m_n * u_m_n + pressure_m, u_m_n * (pressure_m + e_m)};

    // double pressure_p = p[3];
    double e_p = 1/(GAMMA - 1) * pressure_p + 0.5 * rho_p * (u_p_n * u_p_n + u_p_t * u_p_t);
    // Flux flux = {rho * u_n, rho * u_n * u_t, rho * u_n * u_n + pressure, u_n * (pressure + e)};
    Flux flux_p = {rho_p * u_p_n, rho_p * u_p_n * u_p_t, rho_p * u_p_n * u_p_n + pressure_p, u_p_n * (pressure_p + e_p)};

    // ===========================================================
    // A naive implementation: average the fluxes on the both side
    // ===========================================================

    Flux flux = {0.5*(flux_m[0] + flux_p[0]), 0.5*(flux_m[1] + flux_p[1]), 0.5*(flux_m[2] + flux_p[2]), 0.5*(flux_m[3] + flux_p[3])};

    // ===========================================================
    // Roe's average: need another term in the flux
    // ===========================================================
    typedef Eigen::Matrix<double, 4, 4> Matrix4d;
    typedef Eigen::Matrix<double, 4, 1> Vector4d;
    Matrix4d A_abs = Matrix4d::Zero(), Lambda_abs = Matrix4d::Zero();
    auto e_abs = [&](const double &a) { return entropy_abs(a, std::abs(x(1, 0) - x(0, 0)), dt); };
    Matrix4d L = Matrix4d::Zero(), R = Matrix4d::Zero();

    // Private vars here
    {
        double u = u_t_mid, v = u_n_mid;
        double h = h_mid;
        double a = std::sqrt((GAMMA - 1) * (h - 0.5 * (u * u + v * v)));
        Lambda_abs(0, 0) = e_abs(u_n_reconstruct), Lambda_abs(1, 1) = e_abs(u_n_reconstruct), Lambda_abs(2, 2) = e_abs(u_n_reconstruct - a_reconstruct), Lambda_abs(3, 3) = e_abs(u_n_reconstruct + a_reconstruct);

        R << 0, 1, 1, 1,
             1, 0, u, u,
             0, v, v-a, v+a,
             u, 0.5*(v*v-u*u), h-v*a, h+v*a;

        double b2 = (GAMMA - 1) / a / a, b1 = 0.5 * b2 * (u * u + v * v);
        L << -b1*u, 1+b2*u*u, b2*u*v, -b2*u,
             1-b1, b2*u, b2*v, -b2,
             0.5*(b1+v/a), -0.5*b2*u, -0.5*(b2*v+1/a), 0.5*b2,
             0.5*(b1-v/a), -0.5*b2*u, -0.5*(b2*v-1/a), 0.5*b2;
    }

    A_abs = R * Lambda_abs * L;

    // Compute the flux
    Vector4d U_m;
    U_m << m[0], m[0] * u_m_t, m[0] * u_m_n, e_m;
    Vector4d U_p;
    U_p << p[0], p[0] * u_p_t, p[0] * u_p_n, e_p;
    Vector4d Flux_Roe_add = -0.5 * (A_abs * (U_p - U_m));
    flux = {flux[0] + Flux_Roe_add[0], flux[1] + Flux_Roe_add[1], flux[2] + Flux_Roe_add[2], flux[3] + Flux_Roe_add[3]};

    return flux;

}

void Driver::calculate_flux()
{
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 2; i < nx + 2; i++)
    {
        for (int j = 2; j < ny + 2; j++)
        {
            bool at_boundary = j < 4 || j > ny - 1 || i < 4 || i > nx - 1;
            // n average
            {
                double u_n = u(i, j);
                double v_n = v(i, j);

                Vec2 n_bnd = {x(i-1, j-1) - x(i-2, j-1), y(i-1, j-1) - y(i-2, j-1)};
                double len = n_bnd.norm();
                n_bnd = n_bnd.normalize();
                Vec2 n_normal = {-n_bnd.y, n_bnd.x};

                Vec2 u_bnd = {u_n, v_n};
                double u_normal = u_bnd.dot(n_normal);
                double u_tangent = u_bnd.dot(n_bnd);
                
                // double e = 1/(GAMMA - 1) * p_n + 0.5 * rho_n * (u_n * u_n + v_n * v_n);
                // fluxes[i][j].n = {rho_n * u_normal, rho_n * u_normal * u_tangent, rho_n * u_normal * u_normal + p_n, u_normal*(e + p_n)};
                if (at_boundary)
                    fluxes[i][j].n = reconstruct_roe_boundary(get_state(i, j-1), get_state(i, j), get_state(i, j+1), get_state(i, j+2), n_normal);
                else if (use_roe)
                    fluxes[i][j].n = reconstruct_roe(get_state(i, j-1), get_state(i, j), get_state(i, j+1), get_state(i, j+2), n_normal);
                else
                    fluxes[i][j].n = reconstruct_v_l(get_state(i, j-1), get_state(i, j), get_state(i, j+1), get_state(i, j+2), n_normal);
            }

            // e average
            {
                double u_e = u(i, j);
                double v_e = v(i, j);

                Vec2 e_bnd = {x(i-1, j-2) - x(i-1, j-1), y(i-1, j-2) - y(i-1, j-1)};
                double len = e_bnd.norm();
                e_bnd = e_bnd.normalize();
                Vec2 e_normal = {-e_bnd.y, e_bnd.x};

                Vec2 u_bnd = {u_e, v_e};
                double u_normal = u_bnd.dot(e_normal);
                double u_tangent = u_bnd.dot(e_bnd);

                // double e = 1/(GAMMA - 1) * p_e + 0.5 * rho_e * (u_e * u_e + v_e * v_e);
                // fluxes[i][j].e = {rho_e * u_normal, rho_e * u_normal * u_tangent, rho_e * u_normal * u_normal + p_e, u_normal*(e + p_e)};
                if (at_boundary)
                    fluxes[i][j].e = reconstruct_roe_boundary(get_state(i-1, j), get_state(i, j), get_state(i+1, j), get_state(i+2, j), e_normal);
                else if (use_roe)
                    fluxes[i][j].e = reconstruct_roe(get_state(i-1, j), get_state(i, j), get_state(i+1, j), get_state(i+2, j), e_normal);
                else
                    fluxes[i][j].e = reconstruct_v_l(get_state(i-1, j), get_state(i, j), get_state(i+1, j), get_state(i+2, j), e_normal);
            
            }

            // s average
            {
                double u_s = u(i, j);
                double v_s = v(i, j);

                Vec2 s_bnd = {x(i-2, j-2) - x(i-1, j-2), y(i-2, j-2) - y(i-1, j-2)};
                double len = s_bnd.norm();
                s_bnd = s_bnd.normalize();
                Vec2 s_normal = {-s_bnd.y, s_bnd.x};

                Vec2 u_bnd = {u_s, v_s};
                double u_normal = u_bnd.dot(s_normal);
                double u_tangent = u_bnd.dot(s_bnd);

                // double e = 1/(GAMMA - 1) * p_s + 0.5 * rho_s * (u_s * u_s + v_s * v_s);
                // fluxes[i][j].s = {rho_s * u_normal, rho_s * u_normal * u_tangent, rho_s * u_normal * u_normal + p_s, u_normal*(e + p_s)};
                if (at_boundary)
                    fluxes[i][j].s = reconstruct_roe_boundary(get_state(i, j+1), get_state(i, j), get_state(i, j-1), get_state(i, j-2), s_normal);
                else if (use_roe)
                    fluxes[i][j].s = reconstruct_roe(get_state(i, j+1), get_state(i, j), get_state(i, j-1), get_state(i, j-2), s_normal);
                else
                    fluxes[i][j].s = reconstruct_v_l(get_state(i, j+1), get_state(i, j), get_state(i, j-1), get_state(i, j-2), s_normal);

            }

            // w average
            {
                double u_w = u(i, j);
                double v_w = v(i, j);

                Vec2 w_bnd = {x(i-2, j-1) - x(i-2, j-2), y(i-2, j-1) - y(i-2, j-2)};
                double len = w_bnd.norm();
                w_bnd = w_bnd.normalize();
                Vec2 w_normal = {-w_bnd.y, w_bnd.x};

                Vec2 u_bnd = {u_w, v_w};
                double u_normal = u_bnd.dot(w_normal);
                double u_tangent = u_bnd.dot(w_bnd);

                // double e = 1/(GAMMA - 1) * p_w + 0.5 * rho_w * (u_w * u_w + v_w * v_w);
                // fluxes[i][j].w = {rho_w * u_normal, rho_w * u_normal * u_tangent, rho_w * u_normal * u_normal + p_w, u_normal*(e + p_w)};
                if (at_boundary)
                    fluxes[i][j].w = reconstruct_roe_boundary(get_state(i+1, j), get_state(i, j), get_state(i-1, j), get_state(i-2, j), w_normal);
                else if (use_roe)
                    fluxes[i][j].w = reconstruct_roe(get_state(i+1, j), get_state(i, j), get_state(i-1, j), get_state(i-2, j), w_normal);
                else
                    fluxes[i][j].w = reconstruct_v_l(get_state(i+1, j), get_state(i, j), get_state(i-1, j), get_state(i-2, j), w_normal);
            }
        }
    }
}
