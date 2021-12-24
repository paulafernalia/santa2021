// Copyright 2021

#ifndef SANTA2021_H_
#define SANTA2021_H_

#include <iosfwd>
#include <string>
#include <vector>
#include <iomanip>
#include "./gurobi_c++.h"

#define EPS 1e-5
#define GRBVar3D vector<vector<vector<GRBVar>>>
#define GRBVar2D vector<vector<GRBVar>>
#define GRBVar1D vector<GRBVar>
#define intvec vector<int>
#define intvec2D vector<intvec>
#define intvec3D vector<intvec2D>
#define intvec4D vector<intvec3D>

using std::cout;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::cerr;
using std::setprecision;
using std::setw;
using std::setfill;


void
add_constr_one_movie_per_pos(
    int nTeams,
    int nPos,
    int nValues,
    GRBVar3D x,
    GRBModel* pmodel);


void
print_solution(
    int nTeams,
    int nPos,
    int nValues,
    vector<char> names,
    const GRBVar3D& x);


void
print_fractional(
    int nTeams,
    int nPos,
    int nValues,
    vector<char> names,
    const intvec2D& permus,
    const GRBVar3D& x,
    const GRBVar3D& delta,
    const GRBVar2D& gamma);


void
add_vars_x(
    int nTeams,
    int nPos,
    int nValues,
    vector<char> names,
    GRBVar3D* px,
    GRBModel* pmodel);


void
add_vars_length(
    int nPos,
    int nTeams,
    GRBVar* pmaxlength,
    GRBVar1D* pteamlength,
    GRBModel* pmodel);


void
add_vars_delta(
    int nTeams,
    int nPos,
    int nMovies,
    int nPermus,
    GRBVar3D* pdelta,
    GRBModel* pmodel);


void
add_vars_gamma(
    int nTeams,
    int nPermus,
    GRBVar2D* pgamma,
    GRBModel* pmodel);


void
add_constr_max_wildcards(
    int nTeams,
    int nPos,
    int nValues,
    int nWildcards,
    GRBVar3D* px,
    GRBModel* pmodel);


void
add_constr_start_left(
    int nTeams,
    int nPos,
    int nValues,
    int nWildcards,
    GRBVar3D* px,
    GRBModel* pmodel);


void
add_constr_permu_pos(
    int nTeams,
    int nPos,
    int nMovies,
    const intvec2D& permus,
    GRBVar3D* px,
    GRBVar3D* pdelta,
    GRBModel* pmodel);


void
add_constr_permu_team(
    int nTeams,
    int nPos,
    int nMovies,
    const intvec2D& permus,
    GRBVar3D* pdelta,
    GRBVar2D* pgamma,
    GRBModel* pmodel);


void
add_constr_no_permu_if_zero(
    int nTeams,
    int nPos,
    int nMovies,
    const intvec2D& permus,
    GRBVar3D* px,
    GRBVar3D* pdelta,
    GRBModel* pmodel);


void
add_constr_team_symmetry(
    int nTeams,
    int nPos,
    GRBVar3D* px,
    GRBModel* pmodel);


void
add_constr_permu(
    int nTeams,
    const intvec2D& permus,
    GRBVar2D* pgamma,
    GRBModel* pmodel);


int
calculate_nPos(
    int nMovies,
    int nTeams);


vector<char>
generate_names(int nMovies);


void
add_constr_consecutive_movies(
    int nTeams,
    int nPos,
    int nValues,
    GRBVar3D* px,
    GRBModel* pmodel);


void
add_constr_last_movie(
    int nTeams,
    int nPos,
    int nValues,
    GRBVar3D* px,
    GRBVar* pmaxlength,
    GRBVar1D* pteamlength,
    GRBModel* pmodel);


int
factorial(int n) {
    int factorial = 1;

    for (int i = 1; i <= n; ++i) {
        factorial *= i;
    }

    return factorial;
}


intvec2D
generate_all_permus(const vector<int>& values);


vector<int>
generate_all_movies(int nMovies);


void
print_infeasibility(GRBModel* pmodel);

#endif  // SANTA2021_H_
