function v = scale2meanVel(v, meanVel,velMask)
    % v:         [cm/s]
    % meanVel:   [cm/s]
    % v:         Scaled velocity map [cm/s]

    % Ensure meanVel is a vector in the 3rd dimension
    meanVel = reshape(meanVel, 1, 1, []);

    v = v./mean(v(velMask)) .* meanVel;

