// Copyright 2021

#include <iosfwd>
#include <string>
#include <vector>
#include "gurobi_c++.h"

#define EPS 1e-5
#define GRBVar3D vector<vector<vector<GRBVar>>>
#define GRBVar2D vector<vector<GRBVar>>
#define GRBVar1D vector<GRBVar>

using std::cout;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::cerr;


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
add_vars_x(
    int nTeams,
    int nPos,
    int nValues,
    vector<char> names,
    GRBVar3D* px,
    GRBModel* pmodel);


void
add_vars_last_movie(
    int nPos,
    GRBVar* plength,
    GRBModel* pmodel);


void
add_vars_permu_pos(
    int nTeams,
    int nPos,
    int nMovies,
    int nPermus,
    GRBVar3D* pdelta,
    GRBModel* pmodel);


void
add_vars_permu_team(
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
    const vector<vector<int>>& permus,
    GRBVar3D* px,
    GRBVar3D* pdelta,
    GRBModel* pmodel);


void
add_constr_permu_team(
    int nTeams,
    int nPos,
    int nMovies,
    const vector<vector<int>>& permus,
    GRBVar3D* pdelta,
    GRBVar2D* pgamma,
    GRBModel* pmodel);


void
add_constr_permu(
    int nTeams,
    const vector<vector<int>>& permus,
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
    GRBVar* plength,
    GRBModel* pmodel);


int
factorial(int n) {
    int factorial = 1;

    for (int i = 1; i <= n; ++i) {
        factorial *= i;
    }

    return factorial;
}


vector<vector<int>>
generate_all_permus(const vector<int>& values);


vector<int>
generate_all_movies(int nMovies);


void
print_infeasibility(GRBModel* pmodel);
