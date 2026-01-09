# A package for meta-analysis of studies with complex knowledge and unexplained heterogeneity.
# Copyright (C) 2025 Martin Adamčík, e-mail: maths38@gmail.com
# The package requires Python 3.12 or above, with Numpy and Scipy packages, and is OS independent.
# This program comes with ABSOLUTELY NO WARRANTY; for details see the licence.
# This is free software, and you are welcome to redistribute it under certain conditions;
# see the GNU General Public License.

import numpy as np
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
import sys
import re
import itertools

from importlib.resources import files, as_file
import shutil
from pathlib import Path

def copy_sample_input(destination="."):
    src = files("meta_analysis").joinpath("example/input.txt")
    dst = Path(destination) / "input.txt"
    if dst.exists():
        raise FileExistsError(f"{dst} already exists.")
    with as_file(src) as src_path:
        shutil.copy(src_path, dst)

def copy_real_input(destination="."):
    src = files("meta_analysis").joinpath("example/real.txt")
    dst = Path(destination) / "input.txt"
    if dst.exists():
        raise FileExistsError(f"{dst} already exists.")
    with as_file(src) as src_path:
        shutil.copy(src_path, dst)

def main():
    # SETUP - READ AND POSSIBLY CHANGE
    # In any case the procedure converges to a distribution that minimises the weighted sum of
    # cross entropies (Kullback-Leibler divergences) from this distribution to
    # closed convex sets of distributions defined by constraints of the individual studies,
    # where weighting is done proportionally by sample sizes.
    # The resulting set of minimisers in practice often contains a single point.
    # If it does not, then we can tweak which point in this set we approach to by the bias.
    bias = 0  # float(1/1000)
    # Set bias to 0 to get the MaxEnt point if studies are consistent and each contains it.
    # Set bias to float(1/1000) to approach to a point close to the CM infinity point.
    cycles = 0  # Minimal number of iterations of the procedure for testing, 0 is default.
    stopping = 1e-6  # When the distance decrease falls below this we stop the procedure.
    maxiter_const = 1000  # Used in projection as maximum number of iterations.
    # The program will give a warning if maxiter needs to be increased.
    # Maxiter should not be below 500.
    gtol_const = 1e-12  # Used in projection as tolerance for constraint satisfaction.
    # This number determines the maximal level of precision, lower it to 1e-13 if higher precision is needed.
    # This number should not be increased above 1e-11. Use 1e-11 for faster computation and higher stability.
    input_file = "input.txt"  # From where study data are read.
    output_file = "output.txt"  # Where the resulting distribution will be stored.
    log_file = "log.txt"  # Where the iteration process will be logged.
    terminate_for_redundant = True  # The program by default terminates if redundant constraints are detected.
    # The reason is that although often the correct solution is computed, sometimes it is not.
    max_variables = 10  # We cannot have too many variables due to the exponential nature of atomic states.
    # SETUP - READ AND POSSIBLY CHANGE

    # Reading the Input File For Setting Up Variables, Number of Studies, Sample Sizes and Constants
    with open(output_file, "w") as file:  # first we set up an output file
        print(
            'Meta-analysis of studies with complex knowledge and unexplained heterogeneity.\nCopyright (C) 2025 Martin Adamčík, e-mail: maths38@gmail.com\nThe package requires Python 3.12 or above, with Numpy and Scipy packages, and is OS independent.\nThis program comes with ABSOLUTELY NO WARRANTY; for details see the license.\nThis is free software, and you are welcome to redistribute it under certain conditions; see the GNU General Public License.')
        file.write(
            'Meta-analysis of studies with complex knowledge and unexplained heterogeneity.\nCopyright (C) 2025 Martin Adamcik, e-mail: maths38@gmail.com\nThe package requires Python 3.12 or above, with Numpy and Scipy packages, and is OS independent.\nThis program comes with ABSOLUTELY NO WARRANTY; for details see the license.\nThis is free software, and you are welcome to redistribute it under certain conditions; see the GNU General Public License.\n\n')
        file.write('CHECK INPUT:\n')
    with open(log_file, "w") as file:  # second we set up a log file
        file.write(
            'This log file shows how the procedure converged to the output shown in ' + output_file + ' for advanced diagnosis.')
        file.write('\nPlease make sure to resolve any warnings given to you in ' + output_file + '\n')
    try:
        with open(input_file, "r") as file:  # third we try to read the input file
            lines = file.readlines()
        print('Study data were read from', input_file)
        with open(output_file, "a") as file:
            file.write('\nStudy data were read from ' + input_file)
    except FileNotFoundError:
        print('Error: The file', input_file, 'was not found in the directory. The program will terminate.')
        with open(output_file, "a") as file:
            file.write(
                '\nError: The file ' + input_file + ' was not found in the directory. The program will terminate.')
        sys.exit()
    study = 0  # first we assume there are no studies
    variable_pattern = r'\b[A-Za-z_]+\b'  # before the first mention of a study, all variables are listed
    study_pattern = r'study'  # we will ignore the case later, study, STUDY, Study are all accepted
    sample_pattern = r'(\bsample\b|\bsize\b).*?([\d,]+)'
    # sample ... number or size ... number, number can separate thousands with a comma
    variable = []  # here we store our variables
    sample_sizes = []  # here we store sample sizes
    # they need to listed in the file in the same order as the studies
    comment_pattern = r'\#'  # comments in the input file are allowed
    maxiter_pattern = r'\bmaxiter\b.*?([\d,]+)'  # to change maxiter parameter
    stopping_pattern = r'\bstopping\b.*?([\d.]+)'  # to change stopping parameter
    stopping_fraction_pattern = r'\bstopping\b.*?([\d,.]+)\s*/\s*([\d,.]+)'
    stopping_scientific_pattern = r'\bstopping\b.*?1e-([\d]+)'
    gtol_scientific_pattern = r'\bgtol\b.*?1e-([\d]+)'  # to change gtol parameter
    bias_pattern = r'\bbias\b.*?([\d.]+)'  # to change bias parameter
    bias_fraction_pattern = r'\bbias\b.*?([\d,.]+)\s*/\s*([\d,.]+)'
    bias_scientific_pattern = r'\bbias\b.*?1e-([\d]+)'
    for line in lines:
        match = re.search(comment_pattern, line)
        if match:
            parts = re.split(comment_pattern, line)
            line = parts[0]
        match = re.search(study_pattern, line, flags=re.IGNORECASE)
        if match:  # we will know how many studies are there by detecting study keyword
            study += 1
        if study == 0:  # before the first keyword, all words are considered variables
            matches = re.findall(variable_pattern, line)
            variable.extend(matches)
        else:
            match = re.search(sample_pattern, line, flags=re.IGNORECASE)
            if match:  # we capture the sample sizes by dropping commas and converting to int
                try:
                    sample_sizes.append(int(match.group(2).replace(',', '')))
                except ValueError:
                    print("Warning: Could not read a given sample size. The program will terminate.")
                    with open(output_file, "a") as file:
                        file.write('\n\nWarning: Could not read a given sample size. The program will terminate.')
                    sys.exit()
            match = re.search(maxiter_pattern, line, flags=re.IGNORECASE)
            if match:
                try:
                    maxiter_const = int(match.group(1).replace(',', ''))
                except ValueError:
                    maxiter_const = 1000  # default value
                if maxiter_const < 500:
                    maxiter_const = 500  # it must be at least 500
            match = re.search(stopping_pattern, line, flags=re.IGNORECASE)
            if match:
                try:
                    stopping = float(match.group(1))
                except ValueError:
                    stopping = 1e-6
                if stopping > 1e-4:
                    stopping = 1e-4  # it must be at least 1e-4
            match = re.search(stopping_fraction_pattern, line, flags=re.IGNORECASE)
            if match:
                try:
                    stopping = float(match.group(1).replace(',', '')) / float(match.group(2).replace(',', ''))
                except ValueError:
                    stopping = 1e-6
                if stopping > 1e-4:
                    stopping = 1e-4  # it must be at least 1e-4
            match = re.search(stopping_scientific_pattern, line, flags=re.IGNORECASE)
            if match:
                stopping = 1 / (10 ** (int(match.group(1))))
                if stopping > 1e-4:
                    stopping = 1e-4  # it must be at least 1e-4
            match = re.search(bias_pattern, line, flags=re.IGNORECASE)
            if match:
                try:
                    bias = float(match.group(1))
                except ValueError:
                    bias = 0
                if bias > 1:
                    bias = 1
            match = re.search(bias_fraction_pattern, line, flags=re.IGNORECASE)
            if match:
                try:
                    bias = float(match.group(1).replace(',', '')) / float(match.group(2).replace(',', ''))
                except ValueError:
                    bias = 0
                if bias > 1:
                    bias = 1
            match = re.search(bias_scientific_pattern, line, flags=re.IGNORECASE)
            if match:
                bias = 1 / (10 ** (int(match.group(1))))
                if bias > 1:
                    bias = 1
            match = re.search(gtol_scientific_pattern, line, flags=re.IGNORECASE)
            if match:
                gtol_const = 1 / (10 ** (int(match.group(1))))
                if gtol_const < 1e-13:
                    gtol_const = 1e-13
                if gtol_const > 1e-11:
                    gtol_const = 1e-11
    print('The maxiter parameter when computing projections is', maxiter_const,
          'This number must be at least 500, and can be increased if you receive a warning to that end.')
    with open(log_file, "a") as file:
        file.write('\nThe maxiter parameter when computing projections is ' + str(
            maxiter_const) + ' This number must be at least 500, and can be increased if you receive a warning to that end.')
    print('The gtol parameter when computing projections is', gtol_const,
          'You may lower it but not below 1e-13 to increase the maximal achievable level of precision. You may increase it but not above 1e-11 for faster computation and increased stability.')
    with open(log_file, "a") as file:
        file.write('\nThe gtol parameter when computing projections is ' + str(
            gtol_const) + ' You may lower it but not below 1e-13 to increase the maximal achievable level of precision. You may increase it but not above 1e-11 for faster computation and increased stability.\n')
    # Reading the Input File For Setting Up Variables, Number of Studies, Sample Sizes and Constants

    # Setting Up Number of Variables and Number of Studies
    seen = set()  # here we store case-insensitive seen variables
    unique_variable = []  # here we store uniques variables
    for v in variable:
        k = v.casefold()
        if k not in seen:
            seen.add(k)
            unique_variable.append(v)
    variable = unique_variable
    # in the above we removed duplicate variables, where we are case-insensitive in removing
    variables = len(variable)  # number of variables
    if variables == 0:
        print("Warning: No variables inputted. The program will terminate. Make sure to name variables in", input_file,
              "before the first study.")
        with open(output_file, "a") as file:
            file.write(
                '\n\nWarning: No variables inputted. The program will terminate. Make sure to name variables in ' + input_file + ' before the first study.')
        sys.exit()
    if variables > max_variables:
        print("Warning: Too many variables. The program will terminate. We do not support over", max_variables,
              "variables.")
        with open(output_file, "a") as file:
            file.write('\n\nWarning: Too many variables. The program will terminate. We do not support over ' + str(
                max_variables) + ' variables.')
        sys.exit()
    # Check that no keywords are listed as variables
    keywords = {"OR", "NOT", "NEG", "SIZE", "SAMPLE", "MAXITER", "STOPPING", "GTOL", "BIAS"}
    upper_vars = {v.upper() for v in variable}
    if keywords.intersection(upper_vars):
        print(
            "Warning: Keywords used as variables. The program will terminate. Please do not use keywords as variables.")
        with open(output_file, "a") as file:
            file.write(
                '\n\nWarning: Keywords used as variables. The program will terminate. Please do not use keywords as variables.')
        sys.exit()
    studies = study  # number of studies
    if studies != len(sample_sizes):
        print(
            "Warning: The number of sample sizes does not match the number of studies. The program will terminate. Make sure to include one sample size after the word study.")
        with open(output_file, "a") as file:
            file.write(
                '\n\nWarning: The number of sample sizes does not match the number of studies. The program will terminate. Make sure to include one sample size after the word study.')
        sys.exit()
    print('The', len(variable), 'variables observed in', studies, 'studies are', variable)
    with open(output_file, "a") as file:
        file.write(
            '\nThe ' + str(len(variable)) + ' variables observed in ' + str(studies) + ' studies are ' + str(variable))
    # Setting Up Number of Variables and Number of Studies

    # Here we initialise some key parameters using information collected above
    pooled_sample = 0  # to compute pooled sample size
    for j in range(0, studies):
        pooled_sample += sample_sizes[j]
    if pooled_sample == 0:
        print(
            'Warning: The pooled sample size is 0. The program will terminate. At least one sample size must be non-zero.')
        with open(output_file, "a") as file:
            file.write(
                '\n\nWarning: The pooled sample size is 0. The program will terminate. At least one sample size must be non-zero.')
        sys.exit()
    else:
        print('The pooled sample size is', pooled_sample, ' Sample sizes determine study weights.')
        with open(output_file, "a") as file:
            file.write('\nThe pooled sample size is ' + str(pooled_sample) + '. Sample sizes determine study weights.')
    weights = []  # to compute respective study weights
    for j in range(0, studies):
        weights.append(sample_sizes[j] / pooled_sample)
    dimensions = 2 ** variables  # computed number of atoms
    uniform = np.full(dimensions, float(1 / dimensions), dtype=np.float64)  # precompute uniform
    para = uniform.copy()  # main program variable

    # Here we initialise some key parameters using information collected above

    # These procedures transform logical statement to their equivalent sets of atoms
    # conjunction of literals
    # defines a vector of coefficients 'a' to use in linear constraints below
    # variables are positioned 0,1,2,...
    # Y positions of present variables, Y = [], N= [] gives a = [1,1,...], a tautology
    # N positions of variables with negation, Y=[0], N=[0] gives a = [0,0,...], a contradiction
    def conjunction(Y, N):
        a = np.zeros(dimensions, dtype=int)
        for i in range(0, dimensions):
            flag = True
            for j in Y:
                if int(format(i, f'0{variables}b')[j]) != 1:
                    flag = False
            for j in N:
                if int(format(i, f'0{variables}b')[j]) != 0:
                    flag = False
            if flag:
                a[i] = 1
        return a

    # disjunction of conjunctions = normal form
    # defines a vector of coefficients 'a' to use in linear constraints below
    # variables are positioned 0,1,2,...
    # Y[k] positions of present variables in the k-th conjunction
    # N[k] positions of variables with negation in the k-th conjunction
    def normal(Y, N):
        a = np.zeros(dimensions, dtype=int)
        if len(Y) != len(N):
            print(f'Warning: Error in constraints. The program will terminate.')
            with open(output_file, "a") as file:
                file.write('\nWarning: Error in constraints. The program will terminate.')
            sys.exit()
        for k in range(0, len(Y)):
            a = np.logical_or(a, conjunction(Y[k], N[k])).astype(int)
        return a

    # conditional of conjunctions
    def condit_conjunction(YUP, NUP, YCON, NCON, const):
        a = 1 * conjunction(YUP, NUP) - const * conjunction(YCON, NCON)
        return a

    # conditional of normal forms
    def condit_normal(YUP, NUP, YCON, NCON, const):
        a = 1 * normal(YUP, NUP) - const * normal(YCON, NCON)
        return a

    # These procedures transform logical statement to their equivalent sets of atoms

    # Testing the procedures above
    # print('coef', conjunction([0,2],[1]))
    # print('coef', conjunction([2],[]))
    # print('coef', condit_conjunction([0,2],[1], [2], [], 0.2))
    # Testing the procedures above

    # Initialising constraints data, redundant constraints will give warning
    constraints_data = [
        [(np.ones(dimensions), 1.0)]  # sum of all must be one
        for _ in range(studies)
    ]

    # Initialising constraints data, redundant constraints will give warning

    # Each constrain must be on a separate line, and contain variable(s) and a number or a fraction
    # The number has a dot for decimals, but must not have any commas or spaces included
    # The number must be positive and in a decimal from, no scientific notation
    # A fraction can be given using the symbol /, with any or no spaces between it and two numbers
    # For logical AND any symbols, or even no symbols, can be used, except the following
    # For logical OR only these forms can be used: OR, or, Or, can have any number of them
    # For logical NOT only these forms can be used: Not, not, Not, NEG, neg, Neg, !, ~
    # with only any number of spaces before the negated variable name in case of ! ~
    # NOT, NEG, OR are detected as standalone expressions so can appear in variable names
    # If both negated and non-negated forms of a variable appear in conjunction, negated is counted
    # For conditional statement only this symbol can be used: |, and must not have OR included
    # Parentheses are unnecessary and ignored: In order of priority: NOT, AND, OR, CONDITIONAL
    # Variables must not contain numbers, symbols |, /, and are not case-sensitive, but underscores ok
    # Comments are allowed, everything after # is ignored

    def display_constraint(Y, N):
        parts = []
        for y in Y:
            parts.append(variable[y])
        for n in N:
            parts.append(f'NOT {variable[n]}')
        return ' & '.join(parts)

    def reading_conjunction(line):
        y, n = [], []
        for i, var in enumerate(variable):
            # pattern for negation: NOT, NEG, !, or ~ (case-insensitive, optional spaces)
            neg_pattern = rf'(?:\b(?:NOT|NEG)\b|[!~])\s*{re.escape(var)}\b'
            pos_pattern = rf'\b{re.escape(var)}\b'
            # Check for negation
            if re.search(neg_pattern, line, flags=re.IGNORECASE):
                n.append(i)
            # Check for positive occurrence (if not negated)
            elif re.search(pos_pattern, line, flags=re.IGNORECASE):
                y.append(i)
        return y, n

    # Reading the Input File For Constraints
    # Different constraints must be on separate lines
    conditional_pattern = r'\|'  # for conditional statement only this symbol can be used: |
    division_pattern = r'/'  # for division statement only this symbol can be used: \
    or_pattern = r'\bOR\b'
    study = -1
    for line in lines:
        match = re.search(comment_pattern, line)
        if match:
            parts = re.split(comment_pattern, line)
            line = parts[0]
        Y, N = [], []
        C = 1.0
        match = re.search(study_pattern, line, flags=re.IGNORECASE)
        if match:  # we keep track of which study we talk about
            study += 1
            print('Study', study, 'with weight', round(weights[study], 4), 'constraints:')
            with open(output_file, "a") as file:
                file.write(
                    '\n\nStudy ' + str(study) + ' with weight ' + str(round(weights[study], 4)) + ' constraints:\n')
        if study > -1:  # fist, we check for constraints only after first study appeared
            match = re.search(division_pattern, line, flags=re.IGNORECASE)
            if match:
                matching = re.search(r'([\d,.]+)\s*/\s*([\d,.]+)', line)
                if matching:  # we first read two numbers to divide if symbol / is present in between
                    try:
                        C = float(matching.group(1).replace(',', '')) / float(matching.group(2).replace(',', ''))
                    except ValueError:
                        print('Warning: A numerical error in a constraint. The constraint is ignored.')
                        with open(output_file, "a") as file:
                            file.write('\nWarning: A numerical error in a constraint. The constraint is ignored.')
                            continue
            else:
                matching = re.search(r'([\d.]+)', line)
                if matching:  # we first read the single number, which has a dot for decimals
                    try:
                        C = float(matching.group())
                    except ValueError:
                        print('Warning: A numerical error in a constraint. The constraint is ignored.')
                        with open(output_file, "a") as file:
                            file.write('\nWarning: A numerical error in a constraint. The constraint is ignored.')
                            continue
            match = re.search(conditional_pattern, line, flags=re.IGNORECASE)
            if match:  # here we deal with a conditional constraint
                parts = re.split(conditional_pattern, line)
                if len(parts) > 2:
                    print('Warning: Too many conditional symbols | in a constraint. The constraint is ignored.')
                    with open(output_file, "a") as file:
                        file.write(
                            '\nWarning: Too many conditional symbols | in a constraint. The constraint is ignored.')
                else:
                    upper, conditional = parts
                    upper = upper + ' ' + conditional  # upper & conditional is true upper pattern
                    YUP, NUP = reading_conjunction(upper)
                    Y, N = reading_conjunction(conditional)
                    if not (Y == [] and N == []):
                        print('P(', display_constraint(YUP, NUP), ' | ', display_constraint(Y, N), ') = ', C)
                        with open(output_file, "a") as file:
                            file.write('\nP( ' + display_constraint(YUP, NUP) + ' | ' + display_constraint(Y,
                                                                                                           N) + ' ) = ' + str(
                                C))
                        constraints_data[study].append((condit_conjunction(YUP, NUP, Y, N, C), 0.0))
                    else:
                        print('Warning: Empty condition in a constraint. The constraint is ignored.')
                        with open(output_file, "a") as file:
                            file.write('\nWarning: Empty condition in a constraint. The constraint is ignored.')
            else:
                match = re.search(or_pattern, line, flags=re.IGNORECASE)
                if match:  # here we deal with disjunctive constraints
                    parts = re.split(or_pattern, line, flags=re.IGNORECASE)
                    for part in parts:
                        y, n = reading_conjunction(part)
                        Y.append(y)
                        N.append(n)
                    if not (Y == [] and N == []):
                        conj_parts = []
                        for y_sublist, n_sublist in zip(Y, N):
                            conj = display_constraint(y_sublist, n_sublist)
                            conj_parts.append(conj)
                        print('P(', ' OR '.join(conj_parts), ') = ', C)
                        with open(output_file, "a") as file:
                            file.write('\nP( ' + ' OR '.join(conj_parts) + ' ) = ' + str(C))
                        constraints_data[study].append((normal(Y, N), C))
                else:  # finally, the last option is only a single conjunction
                    Y, N = reading_conjunction(line)
                    if not (Y == [] and N == []):
                        print('P(', display_constraint(Y, N), ') = ', C)
                        with open(output_file, "a") as file:
                            file.write('\nP( ' + display_constraint(Y, N) + ' ) = ' + str(C))
                        constraints_data[study].append((conjunction(Y, N), C))
    # Reading the Input File For Constraints

    # Custom constraints data for testing, where variables = 3, studies = 3
    # constraints_data = [
    #    [
    #    (np.ones(dimensions), 1),    # sum of all must be one
    #    (np.array([1, 0, 1, 0, 1, 0, 0, 1]), 0.3),
    #    (np.array([0, 1, 1, 0, 0, 0, 0, 0]), 0.4),  # 0.3-t, 0.4-t, t, 0.3+t
    # ... more (a, c) pairs
    #    ],
    #    [
    #    (np.ones(dimensions), 1),    # sum of all must be one
    #    (np.array([0, 0, 1, 1, 0, 1, 1, 0]), 0.7),
    #    (np.array([0, 1, 1, 0, 0, 1, 0, 0]), 0.2),  # 0.1+q, 0.2-q, q, 0.7-q
    #    (np.array([1, 0, 0, 1, 1, 0, 1, 1]), 0.8),  # 0.1+q, 0.2-q, q, 0.7-q
    # ... more (a, c) pairs
    #    ],
    #    [
    #    (np.ones(dimensions), 1),  # sum of all must be one
    #    (np.array([0, 0, 1, 0, 0, 0, 1, 0]), 0.2),
    #    (np.array([0, 1, 0, 0, 0, 1, 0, 0]), 0.2),
    #    (np.array([0, 0, 1, 1, 1, 1, 0, 0]), 0.1),
    # ... more (a, c) pairs
    #    ]
    # ]
    # Custom constraints data for testing, where variables = 3, studies = 3

    # We need to change constraints data to a different form for optimisation
    constraints = []  # to save all lists of constraints in a required form
    redundancy = False  # redundancy flag
    full_rank = False  # full_rank flag
    for i in range(studies):
        A = np.vstack([a for (a, c) in constraints_data[i]])  # rows of A are the vectors
        c = np.array([c for (a, c) in constraints_data[i]])  # RHS values
        rank_A = np.linalg.matrix_rank(A)
        rank_aug = np.linalg.matrix_rank(np.hstack([A, c.reshape(-1, 1)]))
        if rank_A != rank_aug:
            print(f'Warning: Study {i} has inconsistent constraints. The program will terminate.')
            print(
                f'Make sure to split randomised studies into two: one for the intervention sample and one for the control sample.')
            with open(output_file, "a") as file:
                file.write(f'\n\nWarning: Study {i} has inconsistent constraints. The program will terminate.')
                file.write(
                    f'\n\nMake sure to split randomised studies into two: one for the intervention sample and one for the control sample.')
            sys.exit()
            # inconsistent constraints cannot be processed
        if rank_A < A.shape[0]:
            print(
                f'Warning: Study {i} has redundant constraints. Please remove them to ensure correctness and to increase efficiency.')
            with open(output_file, "a") as file:
                file.write(
                    f'\n\nWarning: Study {i} has redundant constraints. Please remove them to ensure correctness and to increase efficiency.')
            redundancy = True
            # redundant constraints might lead to an incorrectly computed solution, and inefficiency
            if (terminate_for_redundant):
                print(f'The program will terminate.')
                with open(output_file, "a") as file:
                    file.write(f' The program will terminate.')
                sys.exit()  # program will terminate so users would not ignore this
        if rank_A == dimensions:
            full_rank = True
        # print(rank_A, dimensions, full_rank)
        # proj = np.linalg.lstsq(A, c, rcond=None)[0]
        # residual = np.linalg.norm(A @ proj - c)
        # if residual > 1e-8:
        #    print(f'Warning: Study {i} has numerically inconsistent constraints.')
        #    with open(output_file, "a") as file:
        #        file.write(f'\n\nWarning: Study {i} has numerically inconsistent constraints.')
        #    # numerically inconsistent constraints can cause problems as they could not be satisfied using computer arithmetic
        lin_con = LinearConstraint(A, c, c)  # equality: A @ x = c
        constraints.append(lin_con)
    # We need to change constraints data to a different form for optimisation

    # Optimisation (minimisation) procedure for defining the projection of a distribution to
    # a closed convex sets of distribution defined by linear constraints from a given study
    eps = 1e-12  # for numerical stabilisation below, optimiser sometimes uses negative trial values

    # we minimise the Kullback-Leibler divergence (cross entropy) as our objective function
    # def function(x):
    #    return np.sum(x * np.log((x + eps) / (para + eps)))
    def function(x):
        x_clipped = np.clip(x, eps, 1)  # to fix overflow of optimisation procedure
        return np.sum(x_clipped * np.log((x_clipped) / (para)))

    # its gradient
    # def gradient(x):
    #    return np.log((x + eps) / (para + eps)) + 1
    def gradient(x):
        x_clipped = np.clip(x, eps, 1)  # to fix overflow of optimisation procedure
        return np.log((x_clipped) / (para)) + 1

    # its Hessian
    # def hessian(x):
    #    return np.diag(1.0 / x)
    def hessian(x):
        x_clipped = np.clip(x, eps, 1)  # to fix overflow of optimisation procedure
        return np.diag(1.0 / x_clipped)

    maxiter_flag = False  # to check if maxiter was exceeded when computing projections
    def projection(constraints):
        nonlocal maxiter_flag
        # x0 = uniform.copy()  # the starting point is uniform
        x0 = para.copy()  # starting point is current para
        bounds = [(1e-8, 1)] * len(para)  # the lower bound is 1e-8, the upper bound is 1
        if redundancy:  # redundant constraints can be dealt with SVDFactorization, but it can give incorrect solutions
            optioning = {"maxiter": maxiter_const, "gtol": gtol_const, "factorization_method": "SVDFactorization"}
        else:
            optioning = {"maxiter": maxiter_const, "gtol": gtol_const}
        # solving projection of global para into constrained set by minimising the function
        res = minimize(function, x0, jac=gradient, hess=hessian,
                       constraints=constraints, bounds=bounds,
                       method="trust-constr", options=optioning)
        # "gtol": 1e-12 seems to work well, bigger values have problems with increasing distance
        # "maxiter": 500 worked well for less complex meta-analysis, but not for more complex
        if not res.success:
            if "Constraint violation exceeds 'gtol'" in res.message:
                print("Warning: projection failed:", res.message,
                      "Make sure constraints do not force any atom as zero; e.g., use 1 / 100,000 instead of zero. The program will terminate.")  # means a serious problem
                with open(output_file, "a") as file:
                    file.write(
                        "\n\nWarning: projection failed: " + res.message + " Make sure constraints do not force any atom as zero; e.g., use 1 / 100,000 instead of zero. The program will terminate.")
                    sys.exit()
            elif "The maximum number of function evaluations is exceeded." in res.message:
                maxiter_flag = True
            else:
                print("Warning: projection failed:", res.message)
                with open(output_file, "a") as file:
                    file.write(
                        "\n\nWarning: projection failed: " + res.message)
                sys.exit()
        return res.x

    # Optimisation (minimisation) procedure for defining the projection of a distribution to
    # a closed convex sets of distribution defined by linear constraints from a given study

    # The convergence procedure iterates projections with weighted arithmetic pooling
    distance = 10
    previous_distance = 20
    pooled = uniform.copy()  # the initial starting point is uniform
    for i in range(0, cycles):  # first we iterate it for the minimal number of cycles, usually for testing
        para = pooled.copy()  # define new starting point para for projection
        pool = []  # stores all projections for pooling
        for j in range(0, studies):
            pool.append(projection(constraints[j]))
        # pool now contains all projections of para
        pooled = np.full(dimensions, float(0), dtype=np.float64)
        for j in range(0, studies):
            pooled = pooled + pool[j] * weights[j]
        # bias if any
        pooled = pooled * (1 - bias) + bias * uniform
        # the distance must be non-increasing at every iteration
        previous_distance = distance
        distance = 0
        for j in range(0, studies):  # divergences from para to the projections
            distance = distance + function(pool[j]) * weights[j]
        print(f'Iteration {i}:', para, 'Divergence:', distance, 'Difference:', previous_distance - distance)
        with open(log_file, "a") as file:
            file.write('\nIteration ' + str(i) + ': [' + ', '.join(f"{x:.5f}" for x in para) + '] Divergence: ' + str(
                distance) + ' Difference: ' + str(previous_distance - distance) + ' < ' + str(stopping) + ' = stopping')
        if (
                previous_distance - distance) < 0:  # we can get increasing distance if we have reached the programmed projection accuracy level
            print(
                "Warning: Increasing divergence. This could indicate that we may have reached the maximal achievable level of precision. You may lower the default gtol to 1e-13 to increase the maximal achievable level of precision. Type gtol = 1e-13 in the",
                input_file, "after keyword STUDY.")
            with open(log_file, "a") as file:
                file.write(
                    "\n\nWarning: Increasing divergence. This could indicate that we may have reached the maximal achievable level of precision. You may lower the default gtol to 1e-13 to increase the maximal achievable level of precision. Type gtol = 1e-13 in the " + input_file + " after keyword STUDY.")
    counter = cycles
    while (previous_distance - distance) > stopping:  # then we iterate it until we fall below the stopping condition
        para = pooled.copy()  # define new starting point para for projection
        pool = []  # stores all projections for pooling
        for j in range(0, studies):
            pool.append(projection(constraints[j]))
        # pool now contains all projections of para
        pooled = np.full(dimensions, float(0), dtype=np.float64)
        for j in range(0, studies):
            pooled = pooled + pool[j] * weights[j]
        # bias if any
        pooled = pooled * (1 - bias) + bias * uniform
        # the distance must be non-increasing at every iteration
        previous_distance = distance
        distance = 0
        for j in range(0, studies):  # divergences from para to the projections
            distance = distance + function(pool[j]) * weights[j]
        print(f'Iteration {counter}:', para, 'Divergence:', distance, 'Difference:', previous_distance - distance, '<',
              stopping, '= stopping')
        with open(log_file, "a") as file:
            file.write(
                '\nIteration ' + str(counter) + ': [' + ', '.join(f"{x:.5f}" for x in para) + '] Divergence: ' + str(
                    distance) + ' Difference: ' + str(previous_distance - distance) + ' < ' + str(
                    stopping) + ' = stopping')
        counter = counter + 1
        if (
                previous_distance - distance) < 0:  # we can get increasing distance if we have reached the programmed projection accuracy level
            print(
                "Warning: Increasing divergence. This could indicate that we may have reached the maximal achievable level of precision. You may lower the default gtol to 1e-13 to increase the maximal achievable level of precision. Type gtol = 1e-13 in the",
                input_file, "after keyword STUDY.")
            with open(log_file, "a") as file:
                file.write(
                    "\n\nWarning: Increasing divergence. This could indicate that we may have reached the maximal achievable level of precision. You may lower the default gtol to 1e-13 to increase the maximal achievable level of precision. Type gtol = 1e-13 in the " + input_file + " after keyword STUDY.")
            break
    if maxiter_flag:  # if we have exceeded maxiter in any projection
        print(
            "Warning: You may need to increase maxiter to ensure correctness of projections. Type maxiter = 2000 anywhere after keyword STUDY in",
            input_file)
        with open(log_file, "a") as file:
            file.write(
                "\n\nWarning: You may need to increase maxiter to ensure correctness of projections. Type maxiter = 2000 anywhere after keyword STUDY in " + input_file)
    # The convergence procedure iterates projections with weighted arithmetic pooling

    # Individual non-weighted KL-divergences from the resulting distribution to studies are stored here
    values = ['to Study ' + str(j) + ': ' + str(function(pool[j])) for j in range(studies)]
    # Individual non-weighted KL-divergences from the resulting distribution to studies are stored here

    # Final Output
    print('The program finished and its output is written in', output_file)
    with open(output_file, "a") as file:
        file.write('\n\nOUTPUT:\n\n')

    print('The resulting distribution after', counter, 'iterations is', para)
    print(
        'The sum of weighted KL-distances of this distribution to sets of probability functions given by individual studies is:',
        distance)
    print('This sum is what is being minimised. The last improvement to this sum was by ' + str(
        previous_distance - distance) + ' < ' + str(stopping) + ' = stopping.')
    print(
        "Individual non-weighted KL-divergences from the resulting distribution to studies are respectively " + ", ".join(
            values) + ", where smaller values indicate smaller divergences. This could be used to judge the degree to which the individual studies disagree with the resulting distribution, and as a basis of reliability analysis.")

    with open(output_file, "a") as file:
        file.write('The resulting distribution after ' + str(counter) + ' iterations is [' + ', '.join(
            f"{x:.4f}" for x in para) + ']\n')
        file.write(
            'The sum of weighted KL-divergences of this distribution to sets of probability functions given by individual studies is: ' + str(
                distance) + '\n')
        file.write('This sum is what is being minimised. The last improvement to this sum was by ' + str(
            previous_distance - distance) + ' < ' + str(stopping) + ' = stopping\n')

    with open(log_file, "a") as file:
        file.write('\n')

    if studies == 1:
        if full_rank:  # if at least one study's constraint determine a single probability distribution then there is only one optinal solution
            with open(log_file, "a") as file:
                file.write('\nThere is only one optimal solution.\n')
        else:
            if bias == 0:  # with no bias and one study we get maxent
                with open(log_file, "a") as file:
                    file.write(
                        '\nThe resulting distribution is the most entropic point among the distributions consistent with the single given study.\n')
            else:  # or we can go to CM infinity
                with open(log_file, "a") as file:
                    file.write(
                        '\nThe resulting distribution would approach the central mass at infinity distribution among the distributions consistent with the single given study, as the bias set to ' + str(
                            bias) + ' approaches (but is not equal to) 0.\n')
    else:
        if full_rank:  # if at least one study's constraint determine a single probability distribution then there is only one optinal solution
            with open(log_file, "a") as file:
                file.write('\nThere is only one optimal solution.\n')
        else:
            if bias == 0:
                with open(log_file, "a") as file:
                    file.write(
                        '\nIt is unclear which exact point we have approached if there are more optimal solutions. This can be addressed by adding a small bias such as 1/1,000. Write bias = 1/1,000 in the ' + input_file + ' after the keyword STUDY.\n')
            else:
                with open(log_file, "a") as file:
                    file.write(
                        '\nThe resulting distribution would approach the central mass at infinite language distribution among optimal solutions, as the bias set to ' + str(
                            bias) + ' approaches (but is not equal to) 0.\n')

    with open(log_file, "a") as file:
        file.write('\nIndividual non-weighted KL-divergences from the resulting distribution [' + ', '.join(
            f"{x:.4f}" for x in para) + '] to studies are respectively ' + ', '.join(
            values) + ', where smaller values indicate smaller divergences. This could be used to judge the degree to which the individual studies disagree with the resulting distribution, and as a basis of reliability analysis.\n\n')

    with open(output_file, "a") as file:
        file.write('\nAssuming that all available evidence was provided, '
                   'that studies are either homogeneous or it is unknown why they are heterogeneous, '
                   'and the population size is much larger than the pooled sample size, then the following is among the best estimations of probabilities we can make.\n')
        for r in range(1, variables + 1):  # number of variables in the combination
            for combo in itertools.combinations(range(variables), r):
                # For each subset, consider all combinations of negations
                for negation_pattern in itertools.product([False, True], repeat=r):
                    Y = [combo[i] for i in range(r) if not negation_pattern[i]]  # positive
                    N = [combo[i] for i in range(r) if negation_pattern[i]]  # negated
                    labels = " & ".join(
                        ("NOT " if i in N else "") + variable[i] for i in combo
                    )
                    probability = np.sum(conjunction(Y, N) * para)
                    file.write(f"\nP( {labels} ) = {probability:.4f}")

    with open(output_file, "a") as file:
        file.write(
            '\n\nTo find the probability of a disjunction use the above with P( X OR Y ) = P( X ) + P( Y ) - P( X & Y ).')
    with open(output_file, "a") as file:
        file.write(
            '\n\nTo find the probability of a conditional use the above with P( X | Y ) = P( X & Y ) / P( Y ).\n\n')
    # Final Output


if __name__ == "__main__":
    main()
