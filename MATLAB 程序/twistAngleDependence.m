function [ lambda, delta ] = twistAngleDependence( R1, R2, omega, R1_, R2_ , omega_ , mnpq, theta, discardZeroFlag )
%twistAngleDependence implements equations 4 & 7 from the research article
%
% M. Le Ster, T. Märkl & S. Brown, "Moiré Patterns: A Simple Analytical Model"
% 2019, 2D Materials 7, 011005
% https://dx.doi.org/10.1088/2053-1583/ab5470
% 
%	It allows to calculate the wavelength and angle of specific moiré pattern (MP)
%	contributions of two superimposed atomic lattices. They correspond to the
%	difference of two reciprocal space vectors denoted by (m,n) & (p,q) and thus
%	depend on the rotation of two superimposed lattices  R1(_), R2(_), omega(_).
% This function originates from the MPUI program that may serve to explore and
% analyse MPs but it can be used as a stand-alone implementation of the above
% equations. For more information please refer to the MPUI readme file and
% documentation.
%
% USAGE:
% An atomic scale moiré pattern can be described by difference vectors of the
% reciprocal lattices of the two interfering 2-dimensional crystal structures.
% Their real-space structure is defined by R1(_), R2(_), and omega(_), for the
% under- (over-) layer, respectively. The variable theta contains the relative
% twist angles between them. The 4-component vector 'mnpq' contains the indices
% that select one reciprocal lattice vector of each structure. The difference of
% this pair of vectors can be interpreted as one contribution to the MP of the
% two lattices. As a result of varying theta, both the corresponding wavelength
% 'lambda' and moiré pattern fringe angle 'delta' follow a characteristic curve.
%
% NOTE:
% R1(_), R2(_), and omega(_) are all scalars; theta can be a vector of any
% length. Output vectors lambda and delta will have the same length as theta.
% Angles are expected to be in DEGREES and will be output in DEGREES, too.
%
% This function is length-units agnostic. What is chosen as input will be the
% unit for the output of lambda. The user has to take care of any necessary
% conversions!
%
% In special cases, the formula can have a discontinuous point at theta == 0.
% Since the length of the output variables lambda (delta) is supposed to be the
% same as the theta input, we can choose to replace the zero value with NaN to
% avoid the discontinuity. This can be toggled on and off with the
% "discardZeroFlag" (true or false, 1/0).

	if discardZeroFlag
		theta(theta == 0) = NaN;
	end

	kappa1  = mnpq(1) / (R1 * sind(omega ));
	kappa2  = mnpq(2) / (R2 * sind(omega ));
	kappa1_ = mnpq(3) / (R1_* sind(omega_));
	kappa2_ = mnpq(4) / (R2_* sind(omega_));
	
	Delta = kappa1  * kappa2  * cosd(omega) + ...
					kappa1  * kappa1_ * cosd(theta + omega_ - omega) + ... 
					kappa2  * kappa2_ * cosd(theta) + ...
					kappa1_ * kappa2_ * cosd(omega_) - ...
					kappa1  * kappa2_ * cosd(theta - omega) - ...
					kappa2  * kappa1_ * cosd(theta + omega_);

	% fringe period
	lambda = 1./sqrt(kappa1^2 + kappa2^2 + kappa1_^2 + kappa2_^2 - (2 * Delta));
	
	% fringe angle
	delta = atan2d( (kappa1 * cosd(omega) - kappa2 - kappa1_ * cosd(theta + omega_) + kappa2_ * cosd(theta)), ...
					(-kappa1 * sind(omega) + kappa1_ * sind(theta + omega_) - kappa2_ * sind(theta)));
	
	% As before, a delta == 0 discontinuity might occur due to numerical issues.
	% Hence we can suppress this value, too.
	if discardZeroFlag
		delta(delta == 0) = NaN; 
	end
	
	delta = mod(delta, 180)-90; % restrict to +-90degree range

end

