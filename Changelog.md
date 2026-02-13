# Changelog

## Version [1.1.0] - Ambient-pressure-dependent bubble dynamics
### Added

- **Ambient pressure as a tunable input parameter (`p∞`)** in the Gilmore model.
- **Pressure-dependent liquid properties** derived from the Tait equation of state:
  - Ambient mass density `ρ∞(p∞)`
  - Ambient sound speed `c∞(p∞)`
- New example script demonstrating ambient-pressure effects:
   `examples/Ambient_pressure_dependence_example.m`
   `examples/bubble_dyn_acoust_rad_ambient_pressure_example_update_tauL_6ns.m`

------

### Changed

- **Initial conditions updated**:
   `p∞`, `ρ∞`, and `c∞` are now treated as *variables* instead of fixed constants (previously fixed at 1 bar, 20 °C).
- **Extended Gilmore model (`extended_gilmore`) updated** to consistently propagate pressure-dependent:
  - Enthalpy difference between bubble wall and far field
  - Sound speed at the bubble wall
  - Jump-start conditions
- **Kirkwood–Bethe shock calculation (`Kirkwood_Bethe_method`) updated** to account for the modified pressure distribution in the surrounding liquid.
- Pressure inside the bubble remains modeled using an **ideal gas EOS with van der Waals hard core**, but is now indirectly affected by ambient pressure via updated liquid properties.

------

### Physical impact and validation

- In the investigated pressure range **0.1 MPa – 50 MPa**:
  - `ρ∞` increases by ≈ **2.1 %**
  - `c∞` increases by ≈ **6.5 %**
- For a representative bubble wall pressure of **P = 1 GPa**:
  - Enthalpy difference `H` decreases by ≈ **5.6 %**
  - Sound speed at the bubble wall changes only weakly
- These effects lead to **measurable changes in bubble dynamics and shock wave emission**, while remaining physically consistent with previous low-pressure results.

------

### Notes

- Possible pressure dependence of **surface tension** and **viscosity** is currently neglected.
- Results reduce to the previous implementation when `p∞ = 1 bar`, ensuring backward compatibility.



## Version [1.0.0] – Initial Release of LIBDAR

- **LIBDAR** is a Matlab toolbox to calculate laser-induced bubble dynamics, simulate acoustic radiation and track energy partitioning. The acronym LIBDAR stands for "**L**aser-**i**nduced **B**ubble **D**ynamics and **A**coustic **R**adiation".

  Laser-induced cavitation bubbles play a crucial role in laser material processing within liquid environments and in biomedicine and biophotonics, where it enables precise surgery on cells and within transparent tissues. It involves localized energy deposition by laser pulses that results in plasma formation, acoustic radiation, which often evolves into a shock wave, and bubble dynamics. **LIBDAR** serves to provide a comprehensive numerical toolbox to analyze the physical effects of laser-induced bubble dynamics and acoustic radiation.

  The toolbox requires a working base installation of Matlab, with some additional functions from the Signal Processing Toolbox. The toolbox has been tested with Matlab versions 2018b, 2023b and 2024b that are installed on Windows 8, 10 and 11.

  The initial release of **LIBDAR** includes test examples that are able to reproduce the results from the publications: 

  [1] https://doi.org/10.1017/jfm.2022.202;
  
  [2] https://doi.org/10.48550/arXiv.2501.13749; 

  [3] https://doi.org/10.1016/j.ultsonch.2023.106391.

  Further extension of the **LIBDAR** toolbox is on the progress. Please follow the website for the latest release, https://github.com/X-X-Liang/LIBDAR.


  **LIBDAR** is under the copyright of its developers and made available as open-source software under the terms of the GNU General Public License Version 3.0.

