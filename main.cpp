#include <iostream>

#include <cnpy.h>
#include <Eigen/Dense>

#include "Mesh.h"

int main(int argc, char *argv[])
{
    int cnt_input = 100000;
    Driver driver;

    if(argc > 1)
    {
        cnt_input = std::stoi(argv[1]);
    }
    if (argc > 2)
    {
        std::string FluxMethod = argv[2];

        if (FluxMethod == "vl")
        {
            driver.use_roe = false;
        }
        else if(FluxMethod == "roe")
        {
            driver.use_roe = true;
        }
        else
        {
            std::cout << "Invalid flux method. Using Roe's method." << std::endl;
        }
    }
    

    bool save_roe = driver.use_roe;

    driver.apply_initial_conditions();
    driver.apply_boundary_conditions();
    driver.dt = 5e-4;
    bool converge = false;
    int cnt = 0;
    while(!converge && cnt < cnt_input)
    {
        if (cnt < 8000)
        {
            driver.use_roe = true;
        }
        else
        {
            driver.use_roe = save_roe;
            if (cnt % 500 == 0)
            {
                driver.epsilon *= 0.1;
            }
        }
        cnt++;
        // driver.add_viscosity();
        driver.calculate_flux();
        converge = driver.solve(driver.dt);
        // converge = true;
        driver.apply_boundary_conditions();
        if (converge || cnt % 500 == 0)
            std::cout << "Iteration: " << cnt << std::endl;
        // std::cout << driver.rho << std::endl;
        if (cnt < 10)
        {
            converge = false;
        }
    }

    driver.save();

    return 0;
    
}
    