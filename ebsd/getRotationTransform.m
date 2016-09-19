function transMatrix = getRotationTransform(rotationAngle)
% GETROTATIONTRANSFORM   return the euclidean rotation matrix for a given
% angle in radians

transMatrix = [cos(rotationAngle) -sin(rotationAngle); ...
    sin(rotationAngle) cos(rotationAngle)];