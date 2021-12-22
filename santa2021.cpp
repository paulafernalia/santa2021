#include <iostream>
#include <sstream>
#include "santa2021.h"


int
main(int argc,
    char **argv) {

    // Set of movies (#0 is no movie and #1 is the wildcard)
    int nMovies = 3;
    int nTeams = 3;
    int nWildcards = 2;

    int nValues = nMovies + 2;
    vector<char> names = generate_names(nMovies);
    int nPos = calculate_nPos(nMovies, nTeams);
    int nPermus = factorial(nMovies);
    vector<int> movies = generate_all_movies(nMovies);

    vector<vector<int>> permus = generate_all_permus(movies);

    GRBEnv* env = 0;

    // Model
    env = new GRBEnv();
    GRBModel model = GRBModel(*env);
    model.set(GRB_StringAttr_ModelName, "xmasmovies");

    // Define variables
    GRBVar3D x;
    add_vars_x(nTeams, nPos, nValues, names, &x, &model);

    GRBVar3D delta;
    add_vars_permu_pos(nTeams, nPos, nMovies, nPermus, &delta, &model);

    GRBVar2D gamma;
    add_vars_permu_team(nTeams, nPermus, &gamma, &model);

    GRBVar length;
    add_vars_last_movie(nPos, &length, &model);

    // 2) Define duration of the longest schedule
    add_constr_last_movie(nTeams, nPos, nValues, &x, &length, &model);

    // 3) Define constraint one movie in each position
    add_constr_one_movie_per_pos(nTeams, nPos, nValues, x, &model);

    // 4) Start x from left
    add_constr_start_left(nTeams, nPos, nValues, nWildcards, &x,
        &model);

    // 5) No two identical consecutive movies (except when no movie)
    add_constr_consecutive_movies(nTeams, nPos, nValues, &x, &model);

    // 6) Define constraint num wildcards per team
    add_constr_max_wildcards(nTeams, nPos, nValues, nWildcards, &x,
        &model);

    // 4) Find each permutation in each position
    add_constr_permu_pos(nTeams, nPos, nMovies, permus, &x, &delta,
        &model);

    // 5) Find each permutation in each team
    add_constr_permu_team(nTeams, nPos, nMovies, permus, &delta, &gamma,
        &model);

    // 6) Find each permutation anywhere
    add_constr_permu(nTeams, permus, &gamma, &model);


    // Solve problem
    model.write("model.lp");

    model.set(GRB_IntParam_Cuts, 2);
    model.set(GRB_IntParam_Presolve, 2);
    model.optimize();

    // Evaluate result
    int status = model.get(GRB_IntAttr_Status);
    if (status == GRB_INFEASIBLE)
        print_infeasibility(&model);
    else if (status == GRB_OPTIMAL)
        print_solution(nTeams, nPos, nValues, names, x);

    delete env;
    return 0;
}



void
print_infeasibility(GRBModel* pmodel) {
    pmodel->computeIIS();

    // Get constrs
    GRBConstr* constrs = 0;
    constrs = pmodel->getConstrs();

    // Print infeasible
    cout << "\nIIS\n========" << endl;
    for (int c = 0; c < pmodel->get(GRB_IntAttr_NumConstrs); c++) {
        if (constrs[c].get(GRB_IntAttr_IISConstr) == 1)
            cout << constrs[c].get(GRB_StringAttr_ConstrName) << endl;
    }
}


void
add_constr_one_movie_per_pos(
    int nTeams,
    int nPos,
    int nValues,
    GRBVar3D x,
    GRBModel* pmodel) {

    // Define constraint one movie in each position
    for (int g = 0; g < nTeams; g++) {
        // For each position
        for (int t = 0; t < nPos; t++) {
            GRBLinExpr onemovieperpos = 0;
            // For each movie
            for (int m = 0; m < nValues; m++) {
                onemovieperpos += 1 * x[g][t][m];
            }

            // Create constraint name
            ostringstream cname;
            cname << "onemovieperpos(g" << g << ",t" << t << ")";

            // Add constraint
            pmodel->addConstr(onemovieperpos == 1, cname.str());
        }
    }
}


void
print_solution(
    int nTeams,
    int nPos,
    int nValues,
    vector<char> names,
    const GRBVar3D& x) {

    cout << "\n\nSOLUTION" << endl;

    double vsol;
    for (int g = 0; g < nTeams; g++) {
        for (int t = 0; t < nPos; t++) {
            for (int m = 0; m < nValues; m++) {
                // Get solution
                vsol = x[g][t][m].get(GRB_DoubleAttr_X);

                if (vsol > 1 - EPS) {
                    cout << names[m];
                    break;
                }
            }
        }
        cout << endl;
    }
}


void
add_constr_max_wildcards(
    int nTeams,
    int nPos,
    int nValues,
    int nWildcards,
    GRBVar3D* px,
    GRBModel* pmodel) {

    // Add wildcard variable to LHS
    for (int g = 0; g < nTeams; g++) {
        // Create constraint name
        ostringstream cname;
        cname << "maxwildcards(g" << g << ")";

        // Initialise num of wildcards found in team
        GRBLinExpr num_wildcards = 0;

        for (int t = 0; t < nPos; t++)
            num_wildcards += px->at(g)[t][1];

        // Add constraint to model
        pmodel->addConstr(num_wildcards <= nWildcards, cname.str());
    }
}


void
add_constr_start_left(
    int nTeams,
    int nPos,
    int nValues,
    int nWildcards,
    GRBVar3D* px,
    GRBModel* pmodel) {

    for (int g = 0; g < nTeams; g++) {
        for (int t = 0; t < nPos - 1; t++) {
            // Create constraint name
            ostringstream cname;
            cname << "startleft(g" << g << ",t" << t << ")";

            // Add constraint
            pmodel->addConstr(
                px->at(g)[t][0] <= px->at(g)[t + 1][0],
                cname.str());
        }
    }
}


void
add_constr_consecutive_movies(
    int nTeams,
    int nPos,
    int nValues,
    GRBVar3D* px,
    GRBModel* pmodel) {

    for (int g = 0; g < nTeams; g++) {
        for (int t = 0; t < nPos - 1; t++) {
            for (int v = 0; v < nValues; v++) {
                if (v != 0 && v != 1) {
                    // Create constraint name
                    ostringstream cname;
                    cname << "noconsecutive(g" << g << ",t" << t <<
                        ",m" << v << ")";

                    // Add constraint
                    pmodel->addConstr(
                        px->at(g)[t][v] +
                        px->at(g)[t + 1][v] <= 1,
                        cname.str());
                }
            }
        }
    }
}


void
add_constr_permu_pos(
    int nTeams,
    int nPos,
    int nMovies,
    const vector<vector<int>>& permus,
    GRBVar3D* px,
    GRBVar3D* pdelta,
    GRBModel* pmodel) {

    // For each permu
    for (int p = 0; p < permus.size(); p++) {
        // Get this permutation
        const vector<int>& permu = permus[p];

        // For each team
        for (int g = 0; g < nTeams; g++) {
            // For each starting position
            for (int t = 0; t < nPos - nMovies; t++) {
                // For each movie in the permutation
                for (int m = 0; m < permu.size(); m++) {
                    // Create constraint name
                    ostringstream cname;
                    cname << "permutation_position(p" << p << ",g" << g <<
                        ",t" << t << ",m" << m << ")";

                    // Add constraint
                    pmodel->addConstr(
                        pdelta->at(p)[g][t] <=
                        px->at(g)[t + m][permu[m]],
                        cname.str());
                }
            }
        }
    }
}


void
add_constr_permu_team(
    int nTeams,
    int nPos,
    int nMovies,
    const vector<vector<int>>& permus,
    GRBVar3D* pdelta,
    GRBVar2D* pgamma,
    GRBModel* pmodel) {

    // For each permu
    for (int p = 0; p < permus.size(); p++) {
        // For each team
        for (int g = 0; g < nTeams; g++) {
            GRBLinExpr expr = 0;

            // Add all possible positions
            for (int t = 0; t < nPos - nMovies; t++) {
                expr += 1 * pdelta->at(p)[g][t];
            }

            // Create constraint name
            ostringstream cname;
            cname << "permutation_team(p" << p << ",g" << g << ")";

            // Add constraint
            pmodel->addConstr(expr >= pgamma->at(p)[g], cname.str());
        }
    }
}


void
add_constr_permu(
    int nTeams,
    const vector<vector<int>>& permus,
    GRBVar2D* pgamma,
    GRBModel* pmodel) {

    // For each permu
    for (int p = 0; p < permus.size(); p++) {
        const vector<int>& permu = permus[p];

        // If the permu starts with A, B
        if (permu[0] == 2 && permu[1] == 3) {
            // For each team
            for (int g = 0; g < nTeams; g++) {
                // Create constraint name
                ostringstream cname;
                cname << "special_permutation(p" << p << ",g" << g << ")";

                // Add constraint
                pmodel->addConstr(pgamma->at(p)[g] >= 1, cname.str());
            }

        } else {
            GRBLinExpr expr = 0;

            // For each team
            for (int g = 0; g < nTeams; g++) {
                expr += 1 * pgamma->at(p)[g];
            }

            // Create constraint name
            ostringstream cname;
            cname << "regular_permutation(p" << p << ")";

            // Add constraint
            pmodel->addConstr(expr == 1, cname.str());
        }
    }
}

void
add_constr_last_movie(
    int nTeams,
    int nPos,
    int nValues,
    GRBVar3D* px,
    GRBVar* pduration,
    GRBModel* pmodel) {

    // For each team
    for (int g = 0; g < nTeams; g++) {
        GRBLinExpr sum = 0;

        // For each starting position
        for (int t = 0; t < nPos; t++) {
            // For each movie not empty
            for (int v = 0; v < nValues; v++) {
                if (v != 0)
                    sum += px->at(g)[t][v];
            }
        }

        // Create constraint name
        ostringstream cname;
        cname << "length(g" << g << ")";

        // Add constraint
        pmodel->addConstr(*pduration >= sum, cname.str());
    }
}



void
add_vars_x(
    int nTeams,
    int nPos,
    int nValues,
    vector<char> names,
    GRBVar3D* px,
    GRBModel* pmodel) {

    *px = GRBVar3D(nTeams, GRBVar2D());

    // For each team
    for (int g = 0; g < nTeams; g++) {
        px->at(g) = GRBVar2D(nPos, GRBVar1D(nValues));

        // For each position
        for (int t = 0; t < nPos; t++) {
            // For each movie
            for (int m = 0; m < nValues; m++) {
                // Create variable name
                ostringstream vname;
                vname << names[m] << "(g" << g << ",t" << t << ")";

                // Add variable
                px->at(g)[t][m] = pmodel->addVar(0, 1, 0, GRB_BINARY,
                    vname.str());
            }
        }
    }
}


void
add_vars_permu_pos(
    int nTeams,
    int nPos,
    int nMovies,
    int nPermus,
    GRBVar3D* pdelta,
    GRBModel* pmodel) {

    *pdelta = GRBVar3D(nPermus, GRBVar2D());

    // For each permu
    for (int p = 0; p < nPermus; p++) {
        pdelta->at(p) = GRBVar2D(nTeams, GRBVar1D(nPos - nMovies));

        // For each team
        for (int g = 0; g < nTeams; g++) {
            // For each starting position
            for (int t = 0; t < nPos - nMovies; t++) {
                // Create variable name
                ostringstream vname;
                vname << "delta(p" << p << ",g" << g << ",t" << t << ")";

                // Add variable
                pdelta->at(p)[g][t] = pmodel->addVar(0, 1, 0, GRB_BINARY,
                    vname.str());
            }
        }
    }
}


void
add_vars_permu_team(
    int nTeams,
    int nPermus,
    GRBVar2D* pgamma,
    GRBModel* pmodel) {

    *pgamma = GRBVar2D(nPermus, GRBVar1D(nTeams));

    // For each permu
    for (int p = 0; p < nPermus; p++) {
        // For each team
        for (int g = 0; g < nTeams; g++) {
            // Create variable name
            ostringstream vname;
            vname << "gamma(p" << p << ",g" << g << ")";

            // Add variable
            pgamma->at(p)[g] = pmodel->addVar(0, 1, 0, GRB_BINARY,
                vname.str());
        }
    }
}



void
add_vars_last_movie(
    int nPos,
    GRBVar* pduration,
    GRBModel* pmodel) {

    *pduration = pmodel->addVar(0, nPos + 1, 1, GRB_CONTINUOUS, "l");
}


int
calculate_nPos(
    int nMovies,
    int nTeams) {

    float div1, div2;

    if (nMovies <= 2) {
        cerr << "nMovies must be >= 3" << endl;
        exit(3);
    }

    if (nTeams == 1) {
        div1 = 1;
        div2 = 0;
    } else if (nTeams == 2) {
        div1 = 2;
        div2 = 2;
    } else if (nTeams == 3) {
        div1 = nMovies > 3 ? 3 : 2;
        if (nMovies == 3) {
            div2 = 1;
        } else if (nMovies == 4) {
            div2 = 1. / 2.;
        } else {
            div2 = 2. / 3.;
        }
    } else {
        cerr << "nTeams must be 1, 2, 3" << endl;
        exit(3);
    }

    int upper = factorial(nMovies - 1) * (2 * nMovies - 1) / div1 +
        factorial(nMovies - 2) * nMovies * div2;

    cout << "upper = " << upper << endl;

    return upper;
}


vector<char>
generate_names(int nMovies) {
    vector<char> names(nMovies + 2);

    names[0] = '0';
    names[1] = 'W';

    for (int m = 0; m < nMovies; m++)
        names[m + 2] = 65 + m;

    return names;
}


vector<int>
generate_all_movies(int nMovies) {
    vector<int> movies(nMovies);
    for (int m = 0; m < nMovies; m++) {
        movies[m] = m + 2;
    }

    return movies;
}


vector<vector<int>>
generate_all_permus(const vector<int>& values) {
    int nValues = values.size();
    int nPermus = factorial(nValues);

    // Initialise vector to store permutations
    vector<vector<int>> permus(nPermus, vector<int>(nValues));

    int count = 0;

    for (int m = 0; m < nValues; m++) {
        // Eliminate this value from the sequence
        vector<int> subseq(nValues - 1);

        int pos = 0;
        for (int s = 0; s < nValues; s++) {
            if (s != m)
                subseq[pos++] = values[s];
        }

        // Generate subpermutations in a recurrence
        vector<vector<int>> inner_permus = generate_all_permus(subseq);

        // for each subpermutation, append value at the start
        for (int p = 0; p < inner_permus.size(); p++) {
            vector<int> permu = inner_permus[p];
            permu.push_back(values[m]);

            permus[count++] = permu;
        }
    }

    return permus;
}
