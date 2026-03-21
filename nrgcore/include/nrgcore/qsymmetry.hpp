#pragma once
#include <iostream>
#include <string>
#include <vector>

/**
 * @brief Manages quantum number symmetries and conservation laws.
 *
 * qsymmetry is a base class that tracks and manages quantum number symmetries
 * of a system's Hamiltonian. In NRG calculations, symmetries are used to:
 * - Block-diagonalize the Hamiltonian for computational efficiency
 * - Classify states by conserved quantum numbers (e.g., particle number, spin)
 * - Reduce Hilbert space dimension
 *
 * Common symmetries include:
 * - Particle number conservation (charge)
 * - Spin conservation (S_z, S^2)
 * - Orbital angular momentum
 * - Charge parity
 *
 * @example
 * @code
 * qsymmetry sys_sym;
 * sys_sym.add_symmetry("charge");
 * sys_sym.add_symmetry("spin");
 * sys_sym.print_symmetry();  // Lists all symmetries
 * @endcode
 */
class qsymmetry {
  /// @brief Number of symmetries (beyond the base "qsymmetry")
  size_t                   no_of_symmetry{};
  /// @brief List of symmetry identifiers (names of conserved quantum numbers)
  std::vector<std::string> sys_string;

public:
  /**
   * @brief Default constructor initializes base symmetry identifier.
   *
   * Creates an empty symmetry object with only the base "qsymmetry" label.
   * Additional symmetries should be added via add_symmetry().
   */
  qsymmetry() : sys_string({"qsymmetry"}) {}

  /**
   * @brief Add a quantum number symmetry to the system.
   *
   * Registers a new conserved quantum number/symmetry for this system.
   * Examples: "charge", "spin", "particle_number", "parity"
   *
   * @param _id String identifier for the symmetry (should be descriptive)
   *
   * @note Multiple identical symmetries can be added; no uniqueness checking
   * @note Order of addition is preserved for later reference
   */
  void add_symmetry(const std::string &_id) {
    no_of_symmetry++;
    sys_string.push_back(_id);
  }

  /**
   * @brief Get the total number of symmetries (excluding base).
   *
   * @return Number of symmetries that have been added via add_symmetry()
   * @note Does not include the base "qsymmetry" identifier
   */
  [[nodiscard]] size_t get_symmetry_size() const { return no_of_symmetry; }

  /**
   * @brief Print all symmetries to console output.
   *
   * Displays a formatted list of all symmetries including the base identifier.
   * Useful for debugging and verifying symmetry setup.
   *
   * @note Output format: "Symmetry of this model:" followed by indented list
   */
  void                 print_symmetry() const {
    std::cout << "Symmetry of this model:" << std::endl;
    for (const auto &s : sys_string) {
      std::cout << "\t " << s << std::endl;
    }
  }
};
