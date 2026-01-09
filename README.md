# A package for meta-analysis of studies with complex knowledge and unexplained heterogeneity
Copyright (C) 2025 Martin Adamčík, e-mail: maths38@gmail.com  
The package requires Python 3.12 or above, with Numpy and Scipy packages, and is OS independent.
This program comes with ABSOLUTELY NO WARRANTY; for details see the license.
This is free software, and you are welcome to redistribute it under certain conditions; see the GNU General Public License.  


## I. INTRODUCTION

I would like to share with you a tool to combine different studies - to perform their meta-analysis. https://github.com/maths38/meta-analysis_package

In meta-analysis, instead of looking at an individual study, we base our conclusion on several studies thus making it stronger. However, traditional meta-analysis is limited to combining studies that use identical variables. My Python program addresses this limitation by taking into account related information spread across multiple studies weighted by their respective sample sizes.  

If the studies give you the full information on the distribution of their study subjects across the atomic events determined by their sets of variables, then the method is identical to traditional statistical weighted arithmetic mean of these distributions. However, it is more likely that studies are slightly different, and you would need the more sophisticated ability that this program offers. The particular method that this program uses is mathematically proven to be optimal when we combine studies that are more different than we would statistically expect from related studies, but we do not know why. We call this unexplained heterogeneity. On the other hand, if we know why the studies are different, then this (and any other) method is not applicable - in meta-analysis we do not mix apples with oranges!  

Whenever you need to combine studies that involve both complex knowledge and unexplained heterogeneity, this is the optimal extension of the traditional statistical approach. Feel free to share, a GPL license and documentations including the proof of claims are included.  


## II. SETTING UP VARIABLES AND CONSTRAINTS

The method works with probability distributions defined over a finite set of mutually exclusive and exhaustive set of possibilities. This set is defined by the variables that the studies to be combined use to describe observations about their subjects. For example, a study might observe whether patients have cancer, which is the first variable, and suspicious test results, which is the second variable. Lets denote these two variables CANCER and TEST. Then there are four mutually exclusive and exhaustive possibilities for each and every patient: CANCER & TEST, NOT CANCER & TEST, CANCER & NOT TEST, NOT CANCER & NOT TEST. In other words, each patient either has cancer or not, and has suspicious test or not. Ideally, the study also observes exactly how many of their test subjects belong to each of these mutually exclusive categories obtaining a single probability distribution. Nevertheless, especially if several similar studies are combined that do not use the same set of variables, each observation study only gives constraints on the unknown distribution, defined over all variables we work with. 

- For a toy example, say that the first study only reports that  
(a) 7 out of 70 patients had cancer,  
- while the second study reports both that  
(a) 50 out of 100 patients had suspicious result, and  
(b) 20 out of 50 patients who had suspicious test results in fact also had cancer.  

The first study thus says that the probability P of a patient having cancer is P( CANCER ) = 7/70 = 0.1 (or 10%). However, this alone does not determine the exact probability distribution over our mutually exclusive and exhaustive events (also called atoms), only imposes that P (CANCER & TEST ) + P( CANCER & NOT TEST ) = 0.1. Similarly, the second study says that the conditional probability of a patient having cancer given they have suspicious test result is P(CANCER | TEST) = 20/50 = 0.4 and the probability of a patient having suspicious test result is P( TEST ) = 50/100 = 0.5. This, using the formula P( X | Y ) = P( X & Y ) / P( Y ), gives P( CANCER & TEST ) = 0.4*0.5 = 0.2, and from the constraint P( CANCER & TEST ) + P( NOT CANCER & TEST )  = 0.5 we can also infer P( NOT CANCER & TEST ) = 0.3, but we do not know what P( CANCER & NOT TEST ) nor P( CANCER ) precisely are, as they were not measured. However, we have the constraints on what they could be, for example we know that P( CANCER ) is at least 0.2, contradicting the observation from the first study by a level that is statistically unlikely by chance if they selected a sample form the same population. Given we do not have information on what the difference is, we call this unexplained heterogeneity. We can statistically detect heterogeneity using Pearson’s chi-squared test. 

Constraints thus can be and are expected to be inconsistent between different studies, but must be consistent within each and every observation study. It should be theoretically impossible for a study to report contradicting results for a given sample. However, in practice, this commonly happens and the paper Logician Approach to Meta-Analysis with Unexplained Heterogeneity (published in the Journal of Biomedical Informatics 71, https://doi.org/10.1016/j.jbi.2017.05.017) explains how to address the following key considerations. The real meta-analysis described in this paper is encoded in the real example input text file as guidelines.


## III. CONSIDERATIONS NEEDED WHEN SELECTING STUDIES, VARIABLES AND CONSTRAINTS

1. Selecting additional variables increases accuracy by allowing us to compare more apples to apples and oranges to oranges, but it decreases the number of studies that can provide constraints and subjects for these variables thus lowering the overall statistical power. Please also note that the program does not support having over 10 variables.
2. An observation study could report inconsistent constraints when some study subjects leave during the course of the observation. This needs to be addressed by adjusting the reported proportions to make all constraints reported by the study consistent.
3. Experimental studies operate with two samples: control and intervention groups. Such a study must be split into two separate observation studies, each having consistent constraints. Trying to encode an experimental study with one set of constraints would likely lead to an inconsistent set of constraints. It is an inherent weakness of the presented method that it treats an experimental study as two separate observation studies.   
4. Inconsistent constraints can be artificially created by reporting constraints using decimals rather than the original proportions. For example, it is safer to keep 2/9 (reported as 2 out of 9 subjects) rather than using 0.2222. 
5. Constraints must not force any atom to be zero in this method. If a study or definitions of variables require this, use a small constant such as 1 / 100,000 instead of zero. The hard limit below the probability of no atom can go is 1 / 100,000,000.
6. While redundant constraints theoretically should not cause any problems, in this program they do. The program will inform you which study has redundant constraints, and you need to remove them. In the above example, adding P( CANCER & TEST ) = 0.2 would add such a redundant constraint.
7. Reliability analysis can be performed by recomputing the probabilities with studies dropped from the analysis in turns. Dropping a single possibly questionable study should not have a great impact, otherwise the result could be considered unreliable. A simple way to temporarily remove a study is to set its sample size as 0.  


## IV. OPTIMAL SOLUTION

Assuming the population from which we select the samples is much larger than the pooled sample size, and there is unexplained heterogeneity, the following has been proven in On the applicability of the ‘number of possible states’ argument in multi–expert reasoning using a Bayesian approach (published in the Journal of Applied Logic 19, https://doi.org/10.1016/j.jal.2016.10.001). The probability distribution that minimizes the weighted sum of cross entropies (Kullback-Leibler divergences) from this distribution to the closed convex sets of distributions defined by constraints of the individual studies, where weighting is done proportionally to sample sizes, is far more likely to be the distribution found in the population than any other distribution. In the special case when each study determines a single probability distribution, the optimal distribution is the weighted arithmetic mean of the individual studies' distributions, where weighting is again determined proportionally to sample sizes. This agrees with traditional fixed-effect meta-analysis, but extends it to the setting of complex knowledge. On the other hand, this method disagrees with methods of random-effects meta-analysis that manipulate sample weights when heterogeneity is detected.

The shared Python program iteratively produces a sequence of probability distributions that converges to one of the optimal probability distributions defined above. In practice, when many studies with jointly inconsistent constraints are combined, there is likely only one such optimal distribution. However, if all constraints across all studies are jointly consistent, there are many optimal solutions. It is difficult to tell to which exact distribution the process would approach in this situation. However, if we can combine all these consistent constraints into a single study, then the program results in the most entropic distribution that satisfies these constraints, because the iterative procedure starts at the uniform probability distribution. The resulting most entropic distribution is in general considered an optimal representation when we look at one consistent set of constraints that determine a closed convex set of probability distributions. Nevertheless, this is out of the scope of the original justification that focuses on combining jointly inconsistent studies.

If the unknown nature of the distribution to which the procedure approaches in the case of more optimal solutions and multiple studies is unsatisfactory, we can add a small bias: a study claiming the single uniform distribution but with a weight say 1/1,000. If there are more optimal solutions, and regardless of the number of studies, then the procedure would approach the central mass at an infinite language point among these more optimal solutions as the weight approaches (but is not equal to) 0. The proof and explanation can be found in the thesis Collective Reasoning under Uncertainty and Inconsistency (published by the University of Manchester, https://eprints.maths.manchester.ac.uk/id/eprint/2110). In the thesis, we can also find the proofs that this biased procedure, called there LEP or Linear Entropy Process, has the following desirable properties:

1. Language Invariance Principle: This means that we can add variables that we have no constraints for without changing the probabilities restricted to the original language.  
2. Equivalence Principle: The order of studies in which they are arranged does not matter.  
3. Atomic Renaming Principle: If we rename the variables then the resulting probabilities will be renamed accordingly.  
4. Consistency Principle: If all constraints across all studies are jointly consistent, then the procedure approaches a distribution that satisfies all the constraints.  
5. Collegiality Principle: If we add one more study whose constraints are consistent with the previously approached distribution, then this does not change the resulting distribution. In particular, if we add a study with no constraints this will not change the resulting distribution.  
6. Agreement Principle: If two sets of studies result in the same distribution, then adding the studies together would again result in the same distribution.  


## V. UNDERSTANDING OUTPUT

The resulting distribution after several iterations is specified in the output file. Depending on the stopping condition, this distribution may be closer to or further away from one of the optimal solutions defined above. We can increase or decrease the precision (to a certain maximal level) by changing the stopping parameter in the input file or in the source code. The stopping parameter looks at the improvement of the sum of weighted Kullback-Leibler divergences (KL-divergences) of this resulting distribution to sets of probability distributions given by individual studies. During the iteration process, this sum must be non-increasing. An increasing sum could indicate that the optimisation procedure has reached its programmed maximal precision level (given by the gtol parameter), and we cannot improve the accuracy by further changing the stopping condition. Alternatively, this might be a symptom of a problem with the optimisation procedure, and you should resolve any warnings that you may have been given by the program. The last sum of weighted Kullback-Leibler divergences is specified in the output file, and if it is zero (or nearly zero considering the limited precision), then this indicates that the constraints across all studies are jointly consistent and there might be multiple optimal solutions as described above.  

In meta-analysis, the result can be only as good as the input was. The resulting distribution is usable only if all available evidence was provided, that studies are either homogeneous or it is unknown why they are heterogeneous, and population size is much larger than the pooled sample size. If this is the case, you can interpret the resulting distribution by looking at the probabilities of several atomic and non-atomic events included in the output file. The probability of conditionals and disjunctions can be computed from the given probabilities using the following formulas:   

To find the probability of a disjunction use the given probabilities with P( X OR Y ) = P( X ) + P( Y ) - P( X & Y ).  
To find the probability of a conditional use the given probabilities with P( X | Y ) = P( X & Y ) / P( Y ).  

For the toy example discussed above, the program gives the following outputs:  
The resulting distribution is [0.5257, 0.3154, 0.0000, 0.1588]  
P( Cancer ) = 0.1588  
P( NOT Cancer ) = 0.8412  
P( Test ) = 0.4743  
P( NOT Test ) = 0.5257  
P( Cancer & Test ) = 0.1588  
P( Cancer & NOT Test ) = 0.0000  
P( NOT Cancer & Test ) = 0.3154  
P( NOT Cancer & NOT Test ) = 0.5257  
 
Then we can compute that P( CANCER | TEST ) = P( CANCER & TEST | TEST ) = 0.1588 / 0.4743 = 0.3348, which is smaller than what the second study measured (P( CANCER | TEST ) = 20/50 = 0.4). This reflects the information provided by the first study on smaller overall probability of a patient being diagnosed with cancer (P( CANCER ) = 7/70 = 0.1). This way of combining related but different information was not possible to achieve in traditional meta-analysis. The weight of how individual studies contributed to the joint conclusion is fully explainable by their sample sizes - the higher the sample size, the higher the weight.    


## VI. SYNTAX OF INPUT FILE

First list there your propositional variables, and only variables. Variables can use letters and underscores [A-Za-z_].   
STUDY must not be used as a variable even as a part of their name such as STUDYING.  
The following keywords must not be used as variables: SAMPLE, SIZE, MAXITER, STOPPING, GTOL, OR, NOT, NEG, BIAS, but they could be part of their names such as NOTE.  

For the toy example above, we can list variables as follows:

Cancer, Test 

You can separate variables by commas, empty spaces or even put them on separate lines.   

Then the keyword STUDY (not case-sensitive) denotes that new study constraints and the sample size will follow. For each study you must define one sample size. SAMPLE, SIZE and SAMPLE SIZE all can be used as a keyword. For the toy example above, size = 70, sample = 70 or sample size = 70 all can be used. You can add commas to sample size to denote thousands or not. For example, both SAMPLE = 1,250 and SAMPLE = 1250 are fine.  

Now list all constraints for your study on separate lines. They must be consistent. Therefore, make sure to split experimental studies into two: one for the intervention sample and one for the control sample. Additionally, constraints may not force an atomic sentence to be zero. There are three types of constraints: (1) The probability of a conjunction of literals equals a number, (2) the probability of a disjunctive normal form equals a number, and (3) the probability of a conditional equals a number. The probability must be a positive number between 0 and 1. The number has a dot for decimals, but must not have any commas or spaces included. Only the first number included on the line is taken as the probability. The number cannot use the scientific notation, use the decimal form. However, you may want to give probability as a fraction to avoid inconsistencies. If some study subjects leave, this must be ignored to avoid inconsistencies. Use / with any or no spaces between it and two numbers to indicate a fraction. For example, 1/5, 1 / 5, and 1 / 10,000 are all acceptable.  

For the toy example above, we have one simple constraint as follows:

P( CANCER ) =  7/70  

There was no need to use additional symbols, only things that count on a line are variables and the first number. The following alternative forms all could have been used:

P( CANCER ) =  7 / 70  
P(CANCER) =  7 / 70  
P(CANCER) =  0.1  
P CANCER = 7/70  
CANCER = 7/70  
CANCER 7/70    
CANCER 0.1   
  
After we have given the first study's constraints, we need to write the keyword STUDY again to indicate that next study sample size and constraints follow.

STUDY  

sample = 100  

P( TEST ) = 50/100    
P( CANCER | TEST ) = 20/50   

The last constraint of the toy example is a conditional type recognised by key-symbol |. A disjunction cannot appear in conditional constraints, but we can put there any number of conjunctions. For example,

P( CANCER & TEST | TEST ) = 20/50  

is an equivalent constraint to the original. The toy example discussed here is encoded in the toy example input text file as guidelines.


## VII. SYNTAX OF CONSTRAINTS

(1) The probability of a conjunction of literals equals a number

For conjunction (that is logical AND) any symbols, or even no symbols, can be used, except the keywords. We can use any of NOT, NEG, ~, ! for negation. When we use NOT or NEG for negation, we need to add a space unlike with ~ and !. Equivalent examples:

P( NOT CANCER & TEST ) = ...  
P( NEG CANCER & TEST ) = ...  
P( ~CANCER AND TEST ) = ...  
P( !CANCER, TEST ) = ...  

(2) The probability of a disjunctive normal form equals a number

We do not need any parentheses as the priority of logical symbols is as follows: NOT, AND, OR. Therefore the following is equivalent:

P( NOT CANCER & TEST OR CANCER & NOT TEST ) = ...  
P( (NOT CANCER & TEST) OR (CANCER & NOT TEST) ) = ...  
P( ((NOT CANCER) & TEST) OR (CANCER & (NOT TEST)) ) = ...  

Do not write expressions such as P( NOT (CANCER AND TEST) ), NOT will not be recognised. Instead, the disjunctive normal form must be used P( NOT CANCER OR NOT TEST ). There must not be any additional symbols between the negation and the variable except spaces. 

(3) The probability of a conditional equals a number

Equivalent examples:

P( CANCER & TEST | TEST ) = 20/50  
P( CANCER, TEST | TEST ) = 20/50  
P( CANCER TEST | TEST ) = 20/5  

If we had more variables, we could have used logical conjunction also in the condition such as P( Var_A | Var_B & Var_C ). However, conditionals using logical disjunctions (OR) anywhere are not supported.

See the toy and real example input text files for further guidelines.


## VIII. HISTORY
 
This program's origin can be traced back to the Social Entropy Process that George Wilmers proposed for merging individually consistent beliefs of experts that could be jointly inconsistent. The idea was that this would directly extend the Maximum Entropy process that had been already considered the best justified when dealing with individually consistent beliefs of a single expert. The author of the program, as a Ph.D. candidate under George Wilmers's supervision, then tried to use the justification based on counting the number of possible states, which was also used to justify the Maximum Entropy process and was due to Jeff Paris and Alena Vencovská, to justify the Social Entropy Process, but instead found a justification of the first stage of the dual process to the Social Entropy Process, named the Linear Entropy Process. The assumptions used in the proof limit the justification to the case when beliefs of experts are based on observing samples from a population, which size is much larger than the pooled sample size, and the reason as to why the observations are more different than statistically expected is unknown. In a different setting, the original Social Entropy Process due to George Wilmers might be more attractive. The author of the program then used the Linear Entropy Process to combine real medical studies investigating the incidence of cancer with unprovoked venous thromboembolism after suffering from this condition and used the resulting combination as a basis of a real medical advice. This program is being now being publicly shared in hope that others will find it also useful.


## IX. DISCLAIMER

In meta-analysis, the result can be only as good as the input was. The results can be misleading if included studies, variables and constraints were not selected with care. Furthermore, in this program, the resulting distribution is usable only if all available evidence was provided, that studies are either homogeneous or it is unknown why they are heterogeneous, and population size is much larger than the pooled sample size.  


## X. RUNNING THE PROGRAM

1. Install the package using pip:
    ```bash
    pip install meta-analysis
    ```
2. Open a terminal in your working directory (the directory where input and output files will be read and written).
3. Generate a sample toy example `input.txt` in this directory:
   ```bash
   meta-analysis-init
   ``` 
   or the real example `input.txt`:
   ```bash
   meta-analysis-real
   ``` 
4. Edit `input.txt` with the desired study data and parameters.
5. Run:
   ```bash
   meta-analysis
   ``` 
6. Wait until the program finishes.
7. Check `output.txt` in the same directory as `input.txt`.

Alternatively, you may download `input.txt`, and toy and real example text files, from https://github.com/maths38/meta-analysis_package, see the example.

Alternatively, you may run the following Python code from the directory containing `input.txt`.
``` python
from meta_analysis import main
main()
``` 
