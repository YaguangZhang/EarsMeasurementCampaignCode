function [ amps ] = antPatInter(patAz, patEl, azs, els, INTER_METHOD)
%ANTPATINTER Antenna pattern interpolation.
%
% Inputs:
%   - patAz, patEl
%     The antenna patterns, for the Azimuth and Elevation sweeps,
%     respectively; Each of which is a struct containing fields:
%       - azs
%         The azimuth angles in degree from set [0, 360).
%       - els
%         The elevation angles in degree from set [0, 360).
%       - amps
%         The linear amplitudes of the samples.
%       - phases
%         The phases of the samples.
%     All of these fields contains a column vector with each row
%     corresponding to a sweep sample.
%   - azs, els
%     The points that needed to interpolate from angle set [0, 360).
% Output:
%   - amps
%     Interpolated amplitude (linear).
%
% For now, we simply apply a customized linear interpolation.
%
% Yaguang Zhang, Purdue, 10/02/2017

% Methods supported: 'OnLine' and 'WeightedSum'.
if nargin < 5
    INTER_METHOD = 'OnLine';
end

assert(all(size(azs) == size(els)), ...
    'The sizes of inputs azs and els should be the same!');

switch INTER_METHOD
    case 'OnLine'
        %% Linearly Interpolate on the (Azimuth, Elevation) Plane
        
        % Draw a line with slope -1 through the input (az, el) to find the
        % cooresponding crossed points on the x(az) and y(el) axes, i.e.
        %    el = -az + b
        % Note that we only need to worry about the angles ranging from 0
        % to 360 degrees for azs and els, but bs will be ranging from 0 to
        % 720.
        bs = azs+els;
        
        % Fetch the amplitude for the points on the axis via linear
        % interpolation. 
        amp0sAz = interp1([patAz.azs; patAz.azs(2:end)+360], ...
            [patAz.amps; patAz.amps(2:end)], bs);
        amp0sEl = interp1([patEl.els; patEl.els(2:end)+360], ...
            [patEl.amps; patEl.amps(2:end)], bs);
        
        % Find all the 0 b's.
        boolsToSkip = bs(:) == 0;
        indicesToFit = 1:numel(azs);
        indicesToFit = indicesToFit(~boolsToSkip);
        
        % Now linearly interpolate, along the line we drew, between the
        % fetched value, i.e.
        %    (EL0 to P) / (EL0 to AZ0) = az / b 
        %  = (amp - ampEl0) / (ampAz0 ampEl0)
        amps = nan(size(azs));
        amps(indicesToFit) = arrayfun(@(idx) interp1( [0, bs(idx)], ...
            [amp0sEl(idx), amp0sAz(idx)], ...
            azs(idx)), indicesToFit);
        
        % For zero b's, we need to return amp0 for AZ = EL = 0. We will
        % just use the value for zero az from patAz.
        assert(patAz.azs(1) == 0, 'Samples in patAz should start with az=0!');
        amps(boolsToSkip) = patAz.amps(1);
    case 'WeightedSum'
        %% Weighted Sum Metod
        
        % Fetch the corresponding amplitude for (0, el) and (az, 0).
        amp0sAz = interp1(patAz.azs, patAz.amps, azs);
        amp0sEl = interp1(patEl.els, patEl.amps, els);
        % Compute how close the point is to degree 0, and use the result as
        % weights to get the weighted average of the two reference
        % amplitudes, amp2sAz and amp0sEl. Note:
        %     -  360 degree is essentially 0, too;
        %      - We have weights from [0, 1];
        %     -  And the closer the angle is towards 0, the larger the
        %     weight will be;
        %      - The weights for azs are computed via els, and vice versa;
        %      This makes sense because, for example, when az is close to
        %      0, we should depend more on the amplitude from the elevation
        %      sweep.
        wAzs = (cosd(els) + 1)./2;
        wEls = (cosd(azs) + 1)./2;
        amps = (amp0sAz .* wAzs + amp0sEl .* wEls)./(wAzs + wEls);
    otherwise
        error(['Unsupported antenna pattern interpolation method: ', INTER_METHOD, '!'])
        
end

end
% EOF