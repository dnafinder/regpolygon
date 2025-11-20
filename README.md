[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/regpolygon)

📌 Overview
This repository provides the MATLAB function regpolygon, which computes and visualizes properties of regular polygons. A regular polygon is equiangular and equilateral, and this function reports many of its geometric quantities and, for small n, shows how to construct it with ruler-and-compass style steps.

✨ Features
The function regpolygon computes side-based, area-based, and circle-based properties of a regular polygon, including apothem, perimeter, area, circumradius, and more. It reports interior and exterior angles, tests whether the polygon is constructible with classical methods, and displays all results in a MATLAB table. For polygons with 3 to 12 sides, it also provides an animated geometric construction that illustrates how to obtain the regular polygon using ruler and compass. Optional arguments allow you to disable animation and prevent the function from closing existing figures or clearing the Command Window. When requested, regpolygon can return a structured output containing all the computed properties.

🛠 Installation
Download or clone this repository from GitHub:
https://github.com/dnafinder/regpolygon

Add the folder containing regpolygon.m to your MATLAB path using the Add Folder to Path option or the addpath command. No additional toolboxes are required; the function relies only on core MATLAB functionality and basic graphics.

▶️ Usage
You can call regpolygon with positional arguments and optionally control behavior via Name-Value pairs. The simplest use is with the number of sides only, assuming a unit side length. All output is printed in the Command Window and, for certain n, an animated construction is shown in a figure window.

Examples:
regpolygon
regpolygon(5)
regpolygon(8, 2)
regpolygon(7, 1.5, 'Animate', false)
regpolygon(10, 1, 'CloseFigures', false, 'ClearCommandWindow', false)
S = regpolygon(6, 1);  % return a struct of properties

🎛 Inputs
Positional inputs:
n : Number of sides of the regular polygon. Must be an integer greater than or equal to 3. Default is 3.
L : Side length of the polygon. Must be a real scalar greater than 0. Default is 1.

Name-Value options:
'Animate'            : Logical or 0/1 flag. When true, the function uses comet-based animations and arcs to show the construction step by step. When false, only the final polygon, circles, and auxiliary lines are drawn. Default is true.
'CloseFigures'       : Logical or 0/1 flag. When true, all open figures are closed at the beginning of the call. When false, existing figures are preserved. Default is true.
'ClearCommandWindow' : Logical or 0/1 flag. When true, the Command Window is cleared at the beginning to show a clean table and messages. When false, existing output is preserved. Default is true.

📤 Outputs
By default, regpolygon does not return any output variables. Instead, it prints a table in the Command Window listing the main geometric properties of the regular polygon:

Sides
Length
Fixed number f
Apothem
Inscribed circle area
Height
Number of diagonals
Ways of fan triangulations
Perimeter
Area fixed number phi
Area
Circumradius
Circumscribed circle area

It also prints the interior and exterior angles both as rational multiples of pi and in degrees, and indicates whether the polygon is constructible or not.

If you request an output:
S = regpolygon(...)
the function returns a struct S with fields such as S.Sides, S.Length, S.FixedNumber_f, S.Apothem, S.Perimeter, S.Area, S.Circumradius, S.InteriorAngleDeg, S.ExteriorAngleDeg, and S.Constructible. This makes it easy to perform further analyses, create tests, or integrate regpolygon into larger numerical workflows.

🔍 Interpretation
The reported quantities allow you to explore the geometry of regular polygons in a systematic way. The fixed numbers f and phi summarize how area and apothem scale with the side length. The number of diagonals and fan triangulations (based on Catalan numbers) reflect combinatorial properties of the polygon. The constructibility flag is based on classical results in geometry: a regular n-gon is constructible with ruler and compass if and only if n is the product of a power of 2 and distinct Fermat primes.

📝 Notes
For n between 3 and 12, regpolygon also displays an animated geometric construction of the regular polygon. These constructions follow classical or approximate ruler-and-compass techniques. If you set 'Animate' to false, the function draws only the final figure, which can be useful when calling regpolygon repeatedly or when you want a static diagram. For larger n, only the numerical properties are displayed, which is often sufficient for analytical or educational purposes. The constructibility test uses integer factorization of n; for extremely large n this step may become slow.

📚 Citation
If you use this code in scientific, educational, or technical work, please cite it as:

Cardillo G. (2020)
"Regular polygons: facts and drawing".
Available from GitHub:
https://github.com/dnafinder/regpolygon

👤 Author
Author: Giuseppe Cardillo
Email: giuseppe.cardillo.75@gmail.com
GitHub: https://github.com/dnafinder

⚖️ License
This project is distributed under the MIT License. You are free to use, modify, and redistribute the code, provided that the original copyright notice and license text are preserved. The full license terms are provided in the LICENSE file in this GitHub repository.
