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
    intvec movies = generate_all_movies(nMovies);

    intvec2D permus = generate_all_permus(movies);
    intvec3D permu_groups = get_all_permus_by_letter(permus);


    GRBEnv* env = 0;

    // Model
    env = new GRBEnv();
    GRBModel model = GRBModel(*env);
    model.set(GRB_StringAttr_ModelName, "xmasmovies");


    // Define variables
    GRBVar3D x;
    add_vars_x(nTeams, nPos, nValues, names, &x, &model);

    GRBVar3D delta;
    add_vars_delta(nTeams, nPos, nMovies, nPermus, &delta, &model);

    GRBVar2D gamma;
    add_vars_gamma(nTeams, nPermus, &gamma, &model);

    GRBVar maxlength;
    GRBVar1D teamlength;
    add_vars_length(nPos, nTeams, &maxlength, &teamlength, &model);


    // Define constraints

    // 2) Define duration of the longest schedule
    add_constr_last_movie(nTeams, nPos, nValues, &x, &maxlength, &teamlength,
        &model);

    // 3) Define constraint one movie in each position
    add_constr_one_movie_per_pos(nTeams, nPos, nValues, x, &model);

    // 4) Start x from left
    add_constr_start_left(nTeams, nPos, nValues, nWildcards, &x,
        &model);

    // 5) No more than 1 wildcard per subsequence
    add_constr_max_wildcard_sequence(nTeams, nPos, nMovies, &x, &model);

    // 6) Define constraint num wildcards per team
    add_constr_max_wildcards(nTeams, nPos, nValues, nWildcards, &x,
        &model);

    // 7) Find each permutation in each position
    add_constr_permu_pos(nTeams, nPos, nMovies, permus, permu_groups, &x,
        &delta, &model);

    // 8) Find each permutation in each team
    add_constr_permu_team(nTeams, nPos, nMovies, permus, &delta, &gamma,
        &model);

    // 9) 10) Find each permutation anywhere (either globally or within team)
    add_constr_permu(nTeams, permus, &gamma, &model);


    // Define valid inequalities

    // 11) No two identical consecutive movies (except when no movie)
    add_constr_consecutive_movies(nTeams, nPos, nValues, &x, &model);

    // 12) If a zero is found, no sequence covering it counts as a permutation
    add_constr_no_permu_if_zero(nTeams, nPos, nMovies, permus, &x, &delta,
        &model);

    // 13) Break team symmetry
    add_constr_team_symmetry(nTeams, nPos, &x, &model);

    add_constr_two_movies_per_supersequence(nTeams, nPos, nMovies, &x, &model);


    // Solve problem
    model.write("model.lp");

    model.set(GRB_IntParam_Cuts, 2);
    model.set(GRB_IntParam_Presolve, 2);
    model.set(GRB_IntParam_MIPFocus, 1);

    // Solve linear relaxation
    relax_all_vars(&model);
    model.optimize();
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        print_fractional(nTeams, nPos, nValues, names, permus, x, delta, gamma);

    // Solve MIP
    make_all_vars_integer(&model);
    model.optimize();

    // Evaluate MIP result
    int status = model.get(GRB_IntAttr_Status);
    if (status == GRB_INFEASIBLE)
        print_infeasibility(&model);
    else if (status == GRB_OPTIMAL)
        print_solution(nTeams, nPos, nValues, names, x);

    delete env;
    return 0;
}


void
add_constr_max_wildcard_sequence(
    int nTeams,
    int nPos,
    int nMovies,
    GRBVar3D* px,
    GRBModel* pmodel) {

    // For each team and position
    for (int g = 0; g < nTeams; g++) {
        for (int t = 0; t < nPos - nMovies; t++) {
            // Create constraint name
            ostringstream cname;
            cname << "wildcardsequence(g" << g << ",t" << t << ")";

            GRBLinExpr sum = 0;
            for (int m = 0; m < nMovies; m++)
                sum += px->at(g)[t + m][1];

            // Add constraint
            pmodel->addConstr(sum <= 1, cname.str());
        }
    }
}


void
relax_all_vars(
    GRBModel* pmodel) {

    int numVars = pmodel->get(GRB_IntAttr_NumVars);
    GRBVar* vars = pmodel->getVars();

    for (int v = 0; v < numVars; v++)
        vars[v].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
}


void
make_all_vars_integer(
    GRBModel* pmodel) {

    int numVars = pmodel->get(GRB_IntAttr_NumVars);
    GRBVar* vars = pmodel->getVars();

    for (int v = 0; v < numVars; v++)
        vars[v].set(GRB_CharAttr_VType, GRB_INTEGER);
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
print_fractional(
    int nTeams,
    int nPos,
    int nValues,
    vector<char> names,
    const intvec2D& permus,
    const GRBVar3D& x,
    const GRBVar3D& delta,
    const GRBVar2D& gamma) {

    cout << "\n\nX\n=============" << endl;

    double vsol;
    for (int g = 0; g < nTeams; g++) {
        cout << "\nTEAM " << g << ":\n--" << endl;
        for (int v = 0; v < nValues; v++) {
            cout << names[v] << ": ";

            for (int t = 0; t < nPos; t++) {
                // Get solution
                vsol = x[g][t][v].get(GRB_DoubleAttr_X);
                cout << setfill(' ') << setw(5) << setprecision(2) << vsol <<
                    " ";
            }
            cout << endl;
        }
    }


    cout << "\ndelta\n=============" << endl;
    for (int g = 0; g < nTeams; g++) {
        cout << "TEAM " << g << ":\n--" << endl;
        for (int p = 0; p < permus.size(); p++) {
            // cout << "PERMU" << p << ": ";
            for (const auto& letter : permus[p])
                cout << names[letter];
            cout << ": ";

            for (int t = 0; t <= nPos - nValues + 2; t++) {
                // Get solution
                vsol = delta[p][g][t].get(GRB_DoubleAttr_X);
                cout << setfill(' ') << setw(5) << setprecision(2) << vsol <<
                    " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}


void
add_constr_no_permu_if_zero(
    int nTeams,
    int nPos,
    int nMovies,
    const intvec2D& permus,
    GRBVar3D* px,
    GRBVar3D* pdelta,
    GRBModel* pmodel) {

    for (int g = 0; g < nTeams; g++) {
        for (int t = 0; t <= nPos - nMovies; t++) {
            for (int s = 0; s < t + nMovies - 1; s++) {
                GRBLinExpr sum;

                for (int p = 0; p < permus.size(); p++)
                    sum += pdelta->at(p)[g][t];

                sum += px->at(g)[s][0];

                // Create constraint name
                ostringstream cname;
                cname << "nopermuifzero(g" << g << ",t" << t << ",s" << s <<
                    ")";

                // Add constraint
                pmodel->addConstr(sum <= 1, cname.str());
            }
        }
    }
}


void
add_constr_team_symmetry(
    int nTeams,
    int nPos,
    GRBVar3D* px,
    GRBModel* pmodel) {

    for (int t = 0; t < nPos; t++) {
        for (int g = 0; g < nTeams - 1; g++) {
            // Create constraint name
            ostringstream cname;
            cname << "teamsymm(t" << t << ",g" << g << ")";
            // Add constraint
            pmodel->addConstr(
                px->at(g)[t][0] <= px->at(g + 1)[t][0], cname.str());
        }
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
add_constr_two_movies_per_supersequence(
    int nTeams,
    int nPos,
    int nMovies,
    GRBVar3D* px,
    GRBModel* pmodel) {

    // For each team, movie and starting position
    for (int g = 0; g < nTeams; g++) {
        for (int m = 2; m < nMovies + 2; m++) {
            for (int t = 0; t < nPos - 2 * nMovies + 2; t++) {
                // Initialise linear expression
                GRBLinExpr sum = 0;

                // For each intermediate position
                for (int s = t; s < t + 2 * nMovies - 1; s++)
                    sum += px->at(g)[s][m];

                // Create constraint name
                ostringstream cname;
                cname << "max2supersequence(g" << g << ",m" << m << ",t" <<
                    t << ")";

                pmodel->addConstr(sum <= 2, cname.str());
            }
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
            for (int v = 2; v < nValues; v++) {
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


void
add_constr_permu_pos(
    int nTeams,
    int nPos,
    int nMovies,
    const intvec2D& permus,
    const intvec3D& permu_groups,
    GRBVar3D* px,
    GRBVar3D* pdelta,
    GRBModel* pmodel) {

    // For each team
    for (int g = 0; g < nTeams; g++) {
        // For each starting position
        for (int t = 0; t <= nPos - nMovies; t++) {
            // For each relative intermediate position
            for (int s = 0; s < nMovies; s++) {
                // For each movie in the permutation
                for (int m = 2; m < nMovies + 2; m++) {

                    // Define constraint name
                    ostringstream cname;
                    cname << "permutations_position(g" << g << ",m" << m <<
                        ",t" << t << ",s" << s << ")";

                    // Initialise LHS
                    GRBLinExpr sum = 0;

                    // Sum all permutations starting in t with m in position s
                    for (const auto& p : permu_groups[m - 2][s])
                        sum += pdelta->at(p)[g][t];

                    // Add constraint
                    pmodel->addConstr(sum <= px->at(g)[t + s][m],
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
    const intvec2D& permus,
    GRBVar3D* pdelta,
    GRBVar2D* pgamma,
    GRBModel* pmodel) {

    // For each permu
    for (int p = 0; p < permus.size(); p++) {
        // For each team
        for (int g = 0; g < nTeams; g++) {
            GRBLinExpr expr = 0;

            // Add all possible positions
            for (int t = 0; t <= nPos - nMovies; t++)
                expr += pdelta->at(p)[g][t];

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
    const intvec2D& permus,
    GRBVar2D* pgamma,
    GRBModel* pmodel) {

    // For each permu
    for (int p = 0; p < permus.size(); p++) {
        const intvec& permu = permus[p];

        // If the permu starts with A, B
        if (permu[0] == 2 && permu[1] == 3) {
            // For each team
            for (int g = 0; g < nTeams; g++) {
                // Create constraint name
                ostringstream cname;
                cname << "special_permutation(p" << p << ",g" << g << ")";

                // Add constraint
                pmodel->addConstr(pgamma->at(p)[g] == 1, cname.str());
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
            pmodel->addConstr(expr >= 1, cname.str());
        }
    }
}

void
add_constr_last_movie(
    int nTeams,
    int nPos,
    int nValues,
    GRBVar3D* px,
    GRBVar* pmaxlength,
    GRBVar1D* pteamlength,
    GRBModel* pmodel) {

    // For each team
    for (int g = 0; g < nTeams; g++) {
        GRBLinExpr sum = 0;

        // For each starting position
        for (int t = 0; t < nPos; t++) {
            // For each movie not empty
            for (int v = 1; v < nValues; v++) {
                sum += px->at(g)[t][v];
            }
        }

        // Create constraint name
        ostringstream cname_m;
        cname_m << "maxlength(g" << g << ")";

        // Add constraint (max)
        pmodel->addConstr(*pmaxlength >= sum, cname_m.str());

        // Create constraint name
        ostringstream cname_t;
        cname_t << "teamlength(g" << g << ")";

        // Add constraint (max)
        pmodel->addConstr(pteamlength->at(g) == sum, cname_t.str());
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
                px->at(g)[t][m] = pmodel->addVar(0, 1, 0, GRB_CONTINUOUS,
                    vname.str());
            }
        }

        // Add restrictions on empty movies and wilcards
        for (int t = 0; t < nValues - 2; t++) {
            px->at(g)[t][0].set(GRB_DoubleAttr_UB, 0);
            px->at(g)[t][1].set(GRB_DoubleAttr_UB, 0);
        }
    }
}


void
add_vars_delta(
    int nTeams,
    int nPos,
    int nMovies,
    int nPermus,
    GRBVar3D* pdelta,
    GRBModel* pmodel) {

    *pdelta = GRBVar3D(nPermus, GRBVar2D());

    // For each permu
    for (int p = 0; p < nPermus; p++) {
        pdelta->at(p) = GRBVar2D(nTeams, GRBVar1D(nPos - nMovies + 1));

        // For each team
        for (int g = 0; g < nTeams; g++) {
            // For each starting position
            for (int t = 0; t <= nPos - nMovies; t++) {
                // Create variable name
                ostringstream vname;
                vname << "delta(p" << p << ",g" << g << ",t" << t << ")";

                // Add variable
                pdelta->at(p)[g][t] = pmodel->addVar(0, 1, 0, GRB_CONTINUOUS,
                    vname.str());

                pdelta->at(p)[g][t].set(GRB_IntAttr_BranchPriority, 10);
            }
        }
    }
}


void
add_vars_gamma(
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
            pgamma->at(p)[g] = pmodel->addVar(0, 1, 0, GRB_CONTINUOUS,
                vname.str());

            pgamma->at(p)[g].set(GRB_IntAttr_BranchPriority, 100);
        }
    }
}



void
add_vars_length(
    int nPos,
    int nTeams,
    GRBVar* pmaxlength,
    GRBVar1D* pteamlength,
    GRBModel* pmodel) {

    *pmaxlength = pmodel->addVar(0, nPos + 1, 1000, GRB_CONTINUOUS, "ml");
    *pteamlength = vector<GRBVar>(nTeams);

    for (int g = 0; g < nTeams; g++) {
        // Create variable name
        ostringstream vname;
        vname << "teamlength(g" << g << ")";

        pteamlength->at(g) = pmodel->addVar(0, nPos + 1, 1, GRB_CONTINUOUS,
            vname.str());
    }
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


intvec
generate_all_movies(int nMovies) {
    intvec movies(nMovies);
    for (int m = 0; m < nMovies; m++) {
        movies[m] = m + 2;
    }

    return movies;
}


intvec2D
generate_all_permus(const intvec& values) {
    int nValues = values.size();
    int nPermus = factorial(nValues);

    // Initialise vector to store permutations
    intvec2D permus(nPermus, intvec(nValues));

    int count = 0;

    for (int m = 0; m < nValues; m++) {
        // Eliminate this value from the sequence
        intvec subseq(nValues - 1);

        int pos = 0;
        for (int s = 0; s < nValues; s++) {
            if (s != m)
                subseq[pos++] = values[s];
        }

        // Generate subpermutations in a recurrence
        intvec2D inner_permus = generate_all_permus(subseq);

        // for each subpermutation, append value at the start
        for (int p = 0; p < inner_permus.size(); p++) {
            intvec permu = inner_permus[p];
            permu.push_back(values[m]);

            // std::reverse(permu.begin(), permu.end());
            permus[count++] = permu;
        }
    }

    return permus;
}


intvec3D
get_all_permus_by_letter(
    const intvec2D& permus) {

    // Get all unique letters from the first permutation
    intvec unique = permus[0];

    // Initialise vector to store result
    intvec3D permu_groups(unique.size(), intvec2D());

    for (int m = 2; m < unique.size() + 2; m++) {
        permu_groups[m - 2] = intvec2D(unique.size(), intvec());

        for (int t = 0; t < unique.size(); t++) {
            permu_groups[m - 2][t].reserve(unique.size());
        }
    }

    // For each permutation, store where appropriate
    for (int p = 0; p < permus.size(); p++) {
        for (int t = 0; t < unique.size(); t++) {
            int m = permus[p][t];
            permu_groups[m - 2][t].push_back(p);
        }
    }

    return permu_groups;
}
