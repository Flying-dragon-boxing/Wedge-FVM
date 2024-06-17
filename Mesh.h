#include <Eigen/Dense>
#include <string>
#include <vector>

typedef std::array<double, 4> Flux;
typedef std::array<double, 4> State;
typedef std::array<double, 4> Eigenvals;
typedef std::array<Flux, 2> Flux_xy;

struct Node
{
    Flux n, e, s, w;
};

struct Vec2
{
    double x, y;
    double dot(Vec2 v);
    Vec2 operator+(Vec2 v);
    Vec2 operator-(Vec2 v);
    Vec2 operator*(double s);
    Vec2 operator/(double s);
    double norm();
    Vec2 normalize();

};

struct Driver
{
    Driver(Eigen::MatrixXd &x, Eigen::MatrixXd &y);
    Driver();
    ~Driver() = default;

    Eigen::MatrixXd x, y;
    
    int nx, ny;
    double dt = 1e-2;
    double epsilon = 0.08;
    double threshold = 8e-5;

    Eigen::MatrixXd rho, u, v, p;

    std::vector<std::vector<Node>> fluxes;

    // bool use_roe = true;
    std::string scheme = "roe";

    void calculate_flux();
    void apply_boundary_conditions();
    void apply_initial_conditions();
    void add_viscosity();
    double entropy_abs(const double &a, const double &b, const double &c);
    Flux reconstruct_v_l(const State &mm, const State &m, const State &p, const State &pp, Vec2 normal);
    Flux reconstruct_kfvs(const State &mm, const State &m, const State &p, const State &pp, Vec2 normal);
    Flux reconstruct_roe(const State &mm, const State &m, const State &p, const State &pp, Vec2 normal);
    Flux reconstruct_roe_boundary(const State &mm, const State &m, const State &p, const State &pp, Vec2 normal);
    State get_state(int i, int j) {return {rho(i, j), u(i, j), v(i, j), p(i, j)};};
    bool solve(double dt);
    void save();

};
