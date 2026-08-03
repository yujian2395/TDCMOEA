## Introduction
This project contains the MATLAB source code for the Three-Stage Co-Evolutionary Framework Based on Neighborhood Density for Constrained Multi-objective Optimization Problems (**TDCMOEA**). Specifically, the algorithm is designed to solve complex constrained multi-objective optimization problems through a three-stage search strategy, multiple complementary archives, a hybrid dominance relation enhanced by neighborhood density, and diversity enhancement mechanisms integrated throughout the evolutionary process.
## File Descriptions

### 1. **`TDCMOEA.m`**
- **Purpose**: Implements the main logic of the TDCMOEA algorithm.
- **Key Features**:
  - Initializes four complementary archives (`FA`, `FEA`, `DA` and `LPCA`).
  - Groups constraints and progressively activates different constraint groups during the evolutionary process.
  - Executes a three-stage evolutionary process:
    - **First stage**: Performs forward exploration, feasibility exploitation, and diversity enhancement.
    - **Second stage**: Performs feasibility exploitation, diversity enhancement, and layered progressive constraint search.
    - **Third stage**: Focuses on feasibility exploitation and diversity maintenance.
  - Uses the Forward Exploration Archive (`FA`) to explore the unconstrained Pareto front.
  - Uses the Feasibility Exploitation Archive (`FEA`) to search for feasible and convergent solutions.
  - Uses the Diversity Archive (`DA`) to maintain the diversity of the population.
  - Uses the Layered Progressive Constraint Archive (`LPCA`) to progressively handle grouped constraints.

### 2. **`ConstraintGrouping.m`**
- **Purpose**: Groups constraints according to the similarities of their violation patterns.
- **Key Features**:
  - Converts constraint values into binary violation information.
  - Calculates the Jaccard similarity between different constraints.
  - Constructs a distance matrix based on constraint similarities.
  - Divide the constraints into groups
 
### 3. **`DiversityArchive.m`**
- **Purpose**: Maintains an archive of diverse non-dominated solutions.
- **Key Features**:
  - Combines the current archive with newly generated offspring.
  - Calls `DominationCal` to retain non-dominated solutions.
  - Generates uniformly distributed reference vectors.
  - Calculates angular relationships between solutions and reference vectors.
  - Adaptively balances constraint violation and angular diversity according to evolutionary progress and the feasible-solution ratio.
  - Selects representative solutions from different reference-vector neighborhoods.

### 4. **`DominationCal.m`**
- **Purpose**: Calculates the dominance relationships among solutions.
- **Key Features**:
  - Detects objective dominance relationships among solutions.
  - Removes duplicate solutions with identical objective values.
  - Calculates Euclidean distances between solutions in the objective space.
  - Uses local neighborhood density to remove highly similar or crowded solutions.

### 5. **`FeasibilityexploitationArchive.m`**
- **Purpose**: Maintains an archive for feasibility exploitation and convergence improvement.
- **Key Features**:
  - Combines the current archive with newly generated offspring.
  - Calls `DominationCal` to retain non-dominated solutions.
  - Prioritizes feasible solutions according to constraint violation (`CV`).
  - Selects solutions with smaller constraint violations when the number of feasible solutions is insufficient.
  - Normalizes objective values before environmental selection.
  

### 6. **`ForwardexplorationArchive.m`**
- **Purpose**: Maintains an archive for forward exploration in the objective space.
- **Key Features**:
  - Combines the current archive with newly generated offspring.
  - Calls `DominationCal` without constraint-based dominance to retain objective-space non-dominated solutions.
  - Generates uniformly distributed reference vectors.
  - Calculates angular relationships between solutions and reference vectors.
  - Selects representative solutions from different reference-vector neighborhoods.

### 7. **`LayeredProgressiveConstraintArchive.m`**
- **Purpose**: Maintains an archive under progressively activated constraint groups.
- **Key Features**:
  - Evaluates solutions using only the currently active constraints.
  - Filters feasible solutions under the active constraint subset.
  - Calls `DominationCal` to retain non-dominated solutions.
  - Normalizes objective values and generates uniformly distributed reference vectors.
  - Selects representative solutions according to angular relationships.
  - Prioritizes feasible solutions and selects solutions with smaller constraint violations when necessary.

## Instructions for Use

1. Ensure the MATLAB environment is installed and properly configured with the required paths.
2. Place all `.m` files in the same working directory.
3. Use the PlatEMO platform to execute the algorithm. Define the optimization problem using the `Problem` class in PlatEMO.
4. Monitor the optimization results, including the evolution of the population and the final Pareto front.

## Notes

- Input data must conform to the multi-objective optimization problem format, including objective values (`objs`) and constraints (`cons`).
- Algorithm performance depends on parameter settings. The key parameters in the dominance relation function are the distance threshold (`R_eff`) and the neighbor-count difference threshold (`neighborThresh`), which can be fine-tuned in the code.

## Contact
For any questions or suggestions, please contact the author.
