function [ pathLoss ] ...
    = ituSiteSpecificOverRoofTopsNLoS( fInGHz, dInM, ...
    hRInM, wInM, bInM, phiInDegree, h1InM, h2InM, lInM)
%ITUSITESPECIFICOVERROOFTOPSnLOS To compute the path loss for nLoS propagation
%over roof-tops using the site-specific ITU model.
%
% Inputs:
%   - fInGHz
%     The operating frequency in GHz.
%   - dInM
%     3D direct distance between TX and RX.
%   - hRInM
%     The average height of buildings.
%   - wInM
%     The street width.
%   - bInM
%     The average building separation.
%   - phiInDegree
%     The street orientation with respect to the direct path.
%   - h1InM
%     The station 1 antenna height.
%   - h2InM
%     The station 2 antenna height.
%   - lInM
%     The length of the path covered by buildings.
%
% Output:
%   - pathLossInDbGaussian
%     A struct for the resulted path loss, which represents a Gaussian
%     random variable with mean and varience specified by the fields
%     pathLossInDbMean and pathLossInDbVar.
%
% Ref: ITU-R P.1411-9 (06/2017) Annex 1 Section 4.2.2.
%
% Yaguang Zhang, Purdue, 10/16/2017

%% Parameters

F_IN_GHZ_RANG = [0.8, 38];
D_IN_M_RANGE = [10, 5000];

W_IN_M_RANGE = [10, 25];

DELTA_H_1_IN_M_RANGE = [1, 100];
DELTA_H_2_IN_M_RANGE = [4, min(10, hRInM)];

% Make sure the inputs are within the required ranges.
deltaH1InM = h1InM - hRInM;
deltaH2InM = hRInM - h2InM;
if (fInGHz<F_IN_GHZ_RANG(1) || fInGHz>F_IN_GHZ_RANG(2))
    error(['Input fInGHz is out of required range for the ITU model: ', ...
        num2str(F_IN_GHZ_RANG(1)), '~', num2str(F_IN_GHZ_RANG(2))]);
end
if (dInM<D_IN_M_RANGE(1) || dInM>D_IN_M_RANGE(2))
    error(['Input dInM is out of required range for the ITU model: ', ...
        num2str(D_IN_M_RANGE(1)), '~', num2str(D_IN_M_RANGE(2))]);
end
if (wInM<W_IN_M_RANGE(1) || wInM>W_IN_M_RANGE(2))
    error(['Input wInM is out of required range for the ITU model: ', ...
        num2str(W_IN_M_RANGE(1)), '~', num2str(W_IN_M_RANGE(2))]);
end
if (deltaH1InM<DELTA_H_1_IN_M_RANGE(1) ...
        || deltaH1InM>DELTA_H_1_IN_M_RANGE(2))
    error(['Input deltaH1InM is out of required range for the ITU model: ', ...
        num2str(DELTA_H_1_IN_M_RANGE(1)), '~', ...
        num2str(DELTA_H_1_IN_M_RANGE(2))]);
end
if (deltaH2InM<DELTA_H_2_IN_M_RANGE(1) ...
        || deltaH2InM>DELTA_H_2_IN_M_RANGE(2))
    error(['Input deltaH2InM is out of required range for the ITU model: ', ...
        num2str(DELTA_H_2_IN_M_RANGE(1)), '~', ...
        num2str(DELTA_H_2_IN_M_RANGE(2))]);
end

%% Calculation

% Wavelength.
lambdaInM = physconst('LightSpeed')./(fInGHz.*(10.^9));

% We need to find k s.t. d_k <= d <= d_(k+1). For simplicity, just compute
% a bunch of d_k's for comparison.
NUM_KS_TO_EVA = 1000;
ks = 0:(NUM_KS_TO_EVA-1);

% Formula (55).
Aks = wInM.*(h1InM - h2InM).*(2.*ks + 1)./(2.*(hRInM - h2InM));
% Formula (56).
Bks = Aks - ks.*wInM;
% Formula (57).
phiks = atand((Aks./Bks).*tand(phiInDegree));

% Formula (54).
dkps = sqrt( (Aks./sind(phiks)).^(2) + (h1InM - h2InM).^(2) );
% Formula (50).
dks = sqrt( (Bks./sind(phiInDegree)).^(2) + (h1InM - h2InM).^(2) );

% Formula (51).
Ldks = 20.*log10( (4.*pi.*dkps) ./ ( (0.4).^(ks) .* lambdaInM ) );

% Formula (52). Note that the first element in dks is actually d0.
d0 = dks(1);
dksPos = dks(2:end);
dRD = ( 0.25.*dksPos(3) + 0.25.*dksPos(4) ...
    - 0.16.*dksPos(1) - 0.35.*dksPos(2) ) ...
    .*log10(fInGHz) ...
    + 0.25.*dksPos(1) + 0.56.*dksPos(2) ...
    + 0.10.*dksPos(3) + 0.10.*dksPos(4);

% Formula (53). First we need to find k s.t. d_k <= d_RD <= d_(k+1)
%  To be continued...

end
% EOF