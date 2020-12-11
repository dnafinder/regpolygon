# regpolygon
Regular polygon facts and drawing

In Euclidean geometry, a regular polygon is a polygon that is equiangular
(all angles are equal in measure) and equilateral (all sides have the same
length). This function gives all the possible properties of a regulare
polygon. Moreover, for polygons with 3<=n<=12 sides, will be shown an
animation to how contruct it (exactly or approximately) using ruler and
compass. 
 
 Syntax: 	regpolygon(varargin)
    
     Inputs:
           n: sides number, a scalar >=3 (3 by default)
           L: side length, a scalar >0 (1 by default) 
                  
     Outputs:
           Sides
           Length
           Fixed number
           Apotema
           Inscribed circle area
           Height
           number of diagonals
           Perimeter
           Area fixed number (phi)
           Area
           Circumradius
           Circumscribed circle area
           Interior angle 
           Exterior angle 
           Constructible or not 

      Example: 

           Calling on Matlab the function: regpolygon(5)

           Answer is:

                                   Value 
                                  _______
 
     Sides                              5
     Length                             1
     Fixed_Number_(f)             0.68819
     Apotema_(green)              0.68819
     Inscribed_circle_area         1.4879
     Height                        1.5388
     n_of_Diagonals                     5
     Perimeter                          5
     Area_fixed_number_(phi)       1.7205
     Area                          1.7205
     Circumradius_(blue)          0.85065
     Circumscribed_circle_area     2.2733
 
 Interior angle: 3/5*pi	108.00°
 Exterior angle: 2/5*pi	72.00°
 Constructible polygon

           Created by Giuseppe Cardillo
           giuseppe.cardillo.75@gmail.com

 To cite this file, this would be an appropriate format:
 Cardillo G. (2020). Regular polygons: facts and drawing
 https://github.com/dnafinder/regpolygon
