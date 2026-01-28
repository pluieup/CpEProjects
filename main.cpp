#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace nm {

// Utility: clear bad input and flush the current line
static void ClearInputLine()
{
  cin.clear();
  cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

// Utility: robust typed input with a prompt
template <typename T>
static T ReadValue(const string& prompt)
{
  while (true) {
    cout << prompt;
    T value{};
    if (cin >> value) {
      return value;
    }
    cout << "Invalid input. Please try again.\n";
    ClearInputLine();
  }
}

static void WaitForEnter()
{
  cout << "\nPress ENTER to continue...";
  ClearInputLine();
  cin.get();
}

// Format a term like "+3x2" or "-4x3", skipping zero coefficients
static void PrintEquationTermsRow(const vector<double>& row)
{
  const std::size_t n = row.size();
  for (std::size_t j = 0; j < n; ++j) {
    const double a = row[j];
    if (std::abs(a) < 1e-12) {
      continue;
    }

    if (j > 0) {
      if (a >= 0.0) {
        cout << "+";
      }
    }
    cout << a << "x" << (j + 1);
  }
}

// Check if a single equation (row) has diagonal dominance
// i.e., |a_ii| > sum of |a_ij| for all j != i
static bool IsDiagonallyDominant(const vector<double>& row, int diagonalIndex)
{
  if (diagonalIndex < 0 || diagonalIndex >= static_cast<int>(row.size())) {
    return false;
  }

  const double diagonalAbs = std::abs(row[static_cast<std::size_t>(diagonalIndex)]);
  double sumOthers = 0.0;

  for (std::size_t j = 0; j < row.size(); ++j) {
    if (static_cast<int>(j) != diagonalIndex) {
      sumOthers += std::abs(row[j]);
    }
  }

  return diagonalAbs > sumOthers;
}

// Rearrange equations to achieve diagonal dominance using row pivoting
// Returns true if successful, false if diagonal dominance cannot be achieved
static bool RearrangeForDiagonalDominance(vector<vector<double>>& A, vector<double>& b)
{
  const int n = static_cast<int>(A.size());
  vector<int> rowOrder(n);
  for (int i = 0; i < n; ++i) {
    rowOrder[i] = i;
  }

  // Try to assign a row to each diagonal position
  for (int col = 0; col < n; ++col) {
    int bestRow = -1;
    double bestValue = 0.0;

    // Find the row with the largest absolute value in this column
    // among rows not yet assigned
    for (int row = col; row < n; ++row) {
      const int actualRow = rowOrder[row];
      const double absValue = std::abs(A[static_cast<std::size_t>(actualRow)][static_cast<std::size_t>(col)]);
      if (absValue > bestValue) {
        bestValue = absValue;
        bestRow = row;
      }
    }

    if (bestRow == -1 || bestValue < 1e-14) {
      return false; // Cannot find a suitable pivot
    }

    // Swap rows in the order
    if (bestRow != col) {
      std::swap(rowOrder[col], rowOrder[bestRow]);
    }
  }

  // Apply the row permutation
  vector<vector<double>> newA(static_cast<std::size_t>(n));
  vector<double> newB(static_cast<std::size_t>(n));

  for (int i = 0; i < n; ++i) {
    newA[static_cast<std::size_t>(i)] = A[static_cast<std::size_t>(rowOrder[i])];
    newB[static_cast<std::size_t>(i)] = b[static_cast<std::size_t>(rowOrder[i])];
  }

  A = newA;
  b = newB;

  return true;
}

static void Solve_GaussSeidel()
{
  cout << "\n---------------- Gauss-Seidel Method (Topic Set 3) ----------------\n";

  const int n = ReadValue<int>("Enter number of variables (n): ");
  if (n <= 0) {
    cout << "n must be > 0.\n";
    return;
  }

  vector<vector<double>> A(static_cast<std::size_t>(n), vector<double>(static_cast<std::size_t>(n), 0.0));
  vector<double> b(static_cast<std::size_t>(n), 0.0);

  cout << "\nEnter each equation.\n";

  for (int i = 0; i < n; ++i) {
    cout << "Enter equation " << (i + 1) << ": ";
    double value{};
    // First value: constant term b_i
    while (!(cin >> value)) {
      cout << "Invalid input. Please re-enter equation " << (i + 1) << ".\n";
      ClearInputLine();
      cout << "Enter equation " << (i + 1) << ": ";
    }
    b[static_cast<std::size_t>(i)] = value;

    // Next n values: coefficients of x1..xn
    for (int j = 0; j < n; ++j) {
      while (!(cin >> value)) {
        cout << "Invalid coefficient. Please re-enter equation " << (i + 1) << ".\n";
        ClearInputLine();
        cout << "Enter equation " << (i + 1) << ": ";
        j = -1;
        break;
      }
      if (j == -1) {
        // line was re-entered; restart reading full row
        while (!(cin >> value)) {
          cout << "Invalid input. Please re-enter equation " << (i + 1) << ".\n";
          ClearInputLine();
          cout << "Enter equation " << (i + 1) << ": ";
        }
        b[static_cast<std::size_t>(i)] = value;
        j = 0;
      }
      A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = value;
    }
  }

  // Show equations as interpreted
  cout << "\nYour equations:\n";
  for (int i = 0; i < n; ++i) {
    cout << "Equation " << (i + 1) << ": " << b[static_cast<std::size_t>(i)] << " = ";
    PrintEquationTermsRow(A[static_cast<std::size_t>(i)]);
    cout << "\n";
  }

  // Check diagonal dominance before rearrangement
  cout << "\nDiagonal Dominance Check (before rearrangement):\n";
  bool allDominant = true;
  for (int i = 0; i < n; ++i) {
    const bool isDominant = IsDiagonallyDominant(A[static_cast<std::size_t>(i)], i);
    cout << "Equation " << (i + 1) << ": "
         << (isDominant ? "DOMINANT" : "NON-DOMINANT") << "\n";
    if (!isDominant) {
      allDominant = false;
    }
  }

  // Attempt to rearrange for diagonal dominance if needed
  if (!allDominant) {
    cout << "\nAttempting to rearrange equations for diagonal dominance...\n";
    if (RearrangeForDiagonalDominance(A, b)) {
      cout << "Rearrangement successful!\n";
      cout << "\nRearranged equations:\n";
      for (int i = 0; i < n; ++i) {
        cout << "Equation " << (i + 1) << ": " << b[static_cast<std::size_t>(i)] << " = ";
        PrintEquationTermsRow(A[static_cast<std::size_t>(i)]);
        cout << "\n";
      }

      cout << "\nDiagonal Dominance Check (after rearrangement):\n";
      for (int i = 0; i < n; ++i) {
        const bool isDominant = IsDiagonallyDominant(A[static_cast<std::size_t>(i)], i);
        cout << "Equation " << (i + 1) << ": "
             << (isDominant ? "DOMINANT" : "NON-DOMINANT") << "\n";
      }
    } else {
      cout << "Warning: Could not achieve diagonal dominance for all equations.\n";
      cout << "The Gauss-Seidel method may not converge.\n";
    }
  } else {
    cout << "\nAll equations are already diagonally dominant. No rearrangement needed.\n";
  }

  // Show rearranged equations for Gauss-Seidel (solve each for x_i)
  cout << "\nRearranged for Gauss-Seidel iteration:\n";
  for (int i = 0; i < n; ++i) {
    const double aii = A[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
    if (std::abs(aii) < 1e-14) {
      cout << "Warning: diagonal coefficient of equation " << (i + 1)
           << " is zero; method may fail.\n";
    }

    cout << "x" << (i + 1) << " = (" << b[static_cast<std::size_t>(i)];
    for (int j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      const double c = -A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      if (std::abs(c) < 1e-12) {
        continue;
      }
      if (c >= 0.0) {
        cout << "+" << c << "x" << (j + 1);
      } else {
        cout << c << "x" << (j + 1);
      }
    }
    cout << ")/" << aii << "\n";
  }

  cout << "\n--- Iteration Setup ---\n";
  vector<double> x(static_cast<std::size_t>(n), 0.0);
  cout << "Enter initial guess values for x1..x" << n << ":\n";
  for (int i = 0; i < n; ++i) {
    x[static_cast<std::size_t>(i)] =
        ReadValue<double>("x" + std::to_string(i + 1) + " (initial guess): ");
  }

  const double tol = ReadValue<double>("Enter tolerance (e.g., 1e-5): ");
  const int maxIter = ReadValue<int>("Enter max iterations (e.g., 50): ");

  // Calculate precision from tolerance (number of decimal places)
  int precision = 5; // default
  if (tol > 0) {
    precision = -static_cast<int>(std::floor(std::log10(tol)));
  }

  cout << std::fixed << std::setprecision(precision);
  cout << "\nIteration Table:\n";
  cout << "------------------------------------------------------------\n";
  cout << std::left << std::setw(8) << "Iter";
  for (int i = 0; i < n; ++i) {
    cout << std::right << std::setw(12) << ("x" + std::to_string(i + 1));
  }
  cout << "\n";
  cout << "------------------------------------------------------------\n";

  vector<double> xOld = x;
  for (int iter = 1; iter <= maxIter; ++iter) {
    xOld = x;

    for (int i = 0; i < n; ++i) {
      const double aii = A[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
      if (std::abs(aii) < 1e-14) {
        cout << "Diagonal entry of equation " << (i + 1)
             << " became zero. Stopping iterations.\n";
        cout << "------------------------------------------------------------\n";
        return;
      }

      double sum = b[static_cast<std::size_t>(i)];
      for (int j = 0; j < n; ++j) {
        if (j == i) {
          continue;
        }
        sum -= A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
                x[static_cast<std::size_t>(j)];
      }
      x[static_cast<std::size_t>(i)] = sum / aii;
    }

    cout << std::left << std::setw(8) << iter;
    for (int i = 0; i < n; ++i) {
      cout << std::right << std::setw(12) << x[static_cast<std::size_t>(i)];
    }
    cout << "\n";

    double maxDiff = 0.0;
    for (int i = 0; i < n; ++i) {
      maxDiff = std::max(
          maxDiff, std::abs(x[static_cast<std::size_t>(i)] - xOld[static_cast<std::size_t>(i)]));
    }
    if (maxDiff <= tol) {
      break;
    }
  }

  cout << "------------------------------------------------------------\n";
  cout << "End of Gauss-Seidel iterations. Final approximations are in the last row above.\n";
}

} // namespace nm

int main()
{
  using namespace nm;

  cout << "Numerical Methods Computing Machine\n";

  while (true) {
    cout << "\n==================== MAIN MENU ====================\n";
    cout << "Select an option:\n";
    cout << "  [0] Exit\n";
    cout << "  [1] Gauss-Seidel Method (Topic Set 3)\n";

    const int choice = ReadValue<int>("Enter your choice: ");
    switch (choice) {
      case 0:
        cout << "Exiting...\n";
        return 0;
      case 1:
        Solve_GaussSeidel();
        WaitForEnter();
        break;
      default:
        cout << "Invalid choice.\n";
        break;
    }
  }
}

