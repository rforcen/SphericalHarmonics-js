const {
    sqrt
} = Math;

// 3d vector addition
add = (vec1, vec2) => [vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2]];

// 3d scalar multiplication
mult = (c, vec) => [c * vec[0], c * vec[1], c * vec[2]];
// 3d vector subtraction
sub = (vec1, vec2) => [vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2]];

// 3d dot product
dot = (vec1, vec2) =>
    (vec1[0] * vec2[0]) + (vec1[1] * vec2[1]) + (vec1[2] * vec2[2]);

// 3d cross product d1 x d2
cross = (d1, d2) => [(d1[1] * d2[2]) - (d1[2] * d2[1]),
    (d1[2] * d2[0]) - (d1[0] * d2[2]),
    (d1[0] * d2[1]) - (d1[1] * d2[0])
];

// vector norm
mag = vec => sqrt(dot(vec, vec));

// vector magnitude squared
mag2 = vec => dot(vec, vec);

// makes vector unit length
normalize = vec => mult(1 / sqrt(mag2(vec)), vec);

maxv = vec => vec.reduce((pv, cv) => {
    return pv > cv ? pv : cv
})