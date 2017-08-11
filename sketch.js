// physical constants
let g = 9.8;
let m1 = 1;
let m2 = 1;
let l1 = 0.4;
let l2 = 0.6;
let J1 = 1 / 3 * m1 * l1 * l1 / 12;
let J2 = 1 / 3 * m2 * l2 * l2 / 12;

// current variables
var xt = [0, 2 * 3.14 * Math.random(), 2 * 3.14 * Math.random(), 0, 0, 0];
var Et = 0; // energy of pendulum on cart frame (=0 when upright)
var Vt = 0; // Lyapunov function value
var ut = 0; // control input (acceleration of cart)
var dVtdt = 0; // time derivative of the Lyapunov function with linear controller / non-linear model

var current_mode = 'swingup';
var mode_last_change = 0;

// controller related
// linear controller gain
let K = [5.773503e+00, -1.245145e+02, 1.664657e+02, 7.810893e+00, -2.535013e+00, 2.379305e+01, ];
// matrix for Lyapunov function (V(x) = 1/2*x'*S*x)
let S = [
    [1.352886e+00, -4.390771e-01, 4.121077e+00, 9.146507e-01, 3.943857e-01, 7.367118e-01, ],
    [-4.390771e-01, 1.716117e+01, -2.483344e+01, -9.884071e-01, 4.639436e-02, -3.637064e+00, ],
    [4.121077e+00, -2.483344e+01, 5.356611e+01, 4.838637e+00, 1.827593e+00, 8.441639e+00, ],
    [9.146507e-01, -9.884071e-01, 4.838637e+00, 8.617041e-01, 3.812820e-01, 8.374408e-01, ],
    [3.943857e-01, 4.639436e-02, 1.827593e+00, 3.812820e-01, 2.178164e-01, 3.399085e-01, ],
    [7.367118e-01, -3.637064e+00, 8.441639e+00, 8.374408e-01, 3.399085e-01, 1.356479e+00, ],
];

// viewport range
var vxmin = 0,
    vxmax = 0,
    vymin = 0,
    vymax = 0;

// model equation
// dx = f(x) + B(x) u

function f(x) {
    let z = x[0];
    let th1 = x[1];
    let th2 = x[2];
    let dz = x[3];
    let dth1 = x[4];
    let dth2 = x[5];

    var f = [0, 0, 0, 0, 0, 0];

    f[0] = dz;
    f[1] = dth1;
    f[2] = dth2;
    f[4] = ((dth2 * dth2) * l1 * (l2 * l2 * l2) * (m2 * m2) * sin(th1 - th2) * -2.0 + g * l1 * (l2 * l2) * (m2 * m2) * sin(th1) * 4.0 + J2 * g * l1 * m1 * sin(th1) * 8.0 + J2 * g * l1 * m2 * sin(th1) * 1.6E1 - (dth1 * dth1) * (l1 * l1) * (l2 * l2) * (m2 * m2) * cos(th1 - th2) * sin(th1 - th2) * 4.0 - g * l1 * (l2 * l2) * (m2 * m2) * cos(th1 - th2) * sin(th2) * 4.0 + g * l1 * (l2 * l2) * m1 * m2 * sin(th1) * 2.0 - J2 * (dth2 * dth2) * l1 * l2 * m2 * sin(th1 - th2) * 8.0) / (J1 * J2 * 1.6E1 + (l1 * l1) * (l2 * l2) * (m2 * m2) * 4.0 + J2 * (l1 * l1) * m1 * 4.0 + J1 * (l2 * l2) * m2 * 4.0 + J2 * (l1 * l1) * m2 * 1.6E1 - (l1 * l1) * (l2 * l2) * (m2 * m2) * pow(cos(th1 - th2), 2.0) * 4.0 + (l1 * l1) * (l2 * l2) * m1 * m2);
    f[5] = ((dth1 * dth1) * (l1 * l1 * l1) * l2 * (m2 * m2) * sin(th1 - th2) * 8.0 + g * (l1 * l1) * l2 * (m2 * m2) * sin(th2) * 8.0 + J1 * g * l2 * m2 * sin(th2) * 8.0 + (dth1 * dth1) * (l1 * l1 * l1) * l2 * m1 * m2 * sin(th1 - th2) * 2.0 + (dth2 * dth2) * (l1 * l1) * (l2 * l2) * (m2 * m2) * cos(th1 - th2) * sin(th1 - th2) * 4.0 - g * (l1 * l1) * l2 * (m2 * m2) * cos(th1 - th2) * sin(th1) * 8.0 + g * (l1 * l1) * l2 * m1 * m2 * sin(th2) * 2.0 + J1 * (dth1 * dth1) * l1 * l2 * m2 * sin(th1 - th2) * 8.0 - g * (l1 * l1) * l2 * m1 * m2 * cos(th1 - th2) * sin(th1) * 4.0) / (J1 * J2 * 1.6E1 + (l1 * l1) * (l2 * l2) * (m2 * m2) * 4.0 + J2 * (l1 * l1) * m1 * 4.0 + J1 * (l2 * l2) * m2 * 4.0 + J2 * (l1 * l1) * m2 * 1.6E1 - (l1 * l1) * (l2 * l2) * (m2 * m2) * pow(cos(th1 - th2), 2.0) * 4.0 + (l1 * l1) * (l2 * l2) * m1 * m2);

    return f;
}

function B(x) {
    let z = x[0];
    let th1 = x[1];
    let th2 = x[2];
    let dz = x[3];
    let dth1 = x[4];
    let dth2 = x[5];

    var B = [0, 0, 0, 0, 0, 0];
    B[3] = 1.0;
    B[4] = (l1 * (l2 * l2) * (m2 * m2) * cos(th1) * -4.0 - J2 * l1 * m1 * cos(th1) * 8.0 - J2 * l1 * m2 * cos(th1) * 1.6E1 + l1 * (l2 * l2) * (m2 * m2) * cos(th1 - th2) * cos(th2) * 4.0 - l1 * (l2 * l2) * m1 * m2 * cos(th1) * 2.0) / (J1 * J2 * 1.6E1 + (l1 * l1) * (l2 * l2) * (m2 * m2) * 4.0 + J2 * (l1 * l1) * m1 * 4.0 + J1 * (l2 * l2) * m2 * 4.0 + J2 * (l1 * l1) * m2 * 1.6E1 - (l1 * l1) * (l2 * l2) * (m2 * m2) * pow(cos(th1 - th2), 2.0) * 4.0 + (l1 * l1) * (l2 * l2) * m1 * m2);
    B[5] = ((l1 * l1) * l2 * (m2 * m2) * cos(th2) * -8.0 - J1 * l2 * m2 * cos(th2) * 8.0 + (l1 * l1) * l2 * (m2 * m2) * cos(th1 - th2) * cos(th1) * 8.0 - (l1 * l1) * l2 * m1 * m2 * cos(th2) * 2.0 + (l1 * l1) * l2 * m1 * m2 * cos(th1 - th2) * cos(th1) * 4.0) / (J1 * J2 * 1.6E1 + (l1 * l1) * (l2 * l2) * (m2 * m2) * 4.0 + J2 * (l1 * l1) * m1 * 4.0 + J1 * (l2 * l2) * m2 * 4.0 + J2 * (l1 * l1) * m2 * 1.6E1 - (l1 * l1) * (l2 * l2) * (m2 * m2) * pow(cos(th1 - th2), 2.0) * 4.0 + (l1 * l1) * (l2 * l2) * m1 * m2);
    return B;
}

// energy of pendulum on cart frame (=0 when upright)
function E(x) {
    let z = x[0];
    let th1 = x[1];
    let th2 = x[2];
    let dz = x[3];
    let dth1 = x[4];
    let dth2 = x[5];

    let Ex = J1 * (dth1 * dth1) * (1.0 / 2.0) + J2 * (dth2 * dth2) * (1.0 / 2.0) + (dth1 * dth1) * (l1 * l1) * m1 * (1.0 / 8.0) + (dth1 * dth1) * (l1 * l1) * m2 * (1.0 / 2.0) + (dth2 * dth2) * (l2 * l2) * m2 * (1.0 / 8.0) - g * l1 * m1 * (1.0 / 2.0) - g * l1 * m2 - g * l2 * m2 * (1.0 / 2.0) + g * l1 * m1 * cos(th1) * (1.0 / 2.0) + g * l1 * m2 * cos(th1) + g * l2 * m2 * cos(th2) * (1.0 / 2.0) + dth1 * dth2 * l1 * l2 * m2 * cos(th1 - th2) * (1.0 / 2.0);
    return Ex;
}

// time derivative of energy / control input
function dEc(x) {
    let z = x[0];
    let th1 = x[1];
    let th2 = x[2];
    let dz = x[3];
    let dth1 = x[4];
    let dth2 = x[5];

    dExc = dth1 * l1 * m1 * cos(th1) * (-1.0 / 2.0) - dth1 * l1 * m2 * cos(th1) - dth2 * l2 * m2 * cos(th2) * (1.0 / 2.0);
    return dExc;
}



function dVdt(x) {
    let z = x[0];
    let th1 = x[1];
    let th2 = x[2];
    let dz = x[3];
    let dth1 = x[4];
    let dth2 = x[5];

    let dVxdt = -(S[3][4] * dth1 + S[3][5] * dth2 + S[3][3] * dz + S[1][3] * th1 + S[2][3] * th2 + S[0][3] * z) * (K[4] * dth1 + K[5] * dth2 + K[3] * dz + K[1] * th1 + K[2] * th2 + K[0] * z) + dth1 * (S[1][4] * dth1 + S[1][5] * dth2 + S[1][3] * dz + S[1][1] * th1 + S[1][2] * th2 + S[0][1] * z) + dth2 * (S[2][4] * dth1 + S[2][5] * dth2 + S[2][3] * dz + S[1][2] * th1 + S[2][2] * th2 + S[0][2] * z) + dz * (S[0][4] * dth1 + S[0][5] * dth2 + S[0][3] * dz + S[0][1] * th1 + S[0][2] * th2 + S[0][0] * z) + (l1 * (S[4][4] * dth1 + S[4][5] * dth2 + S[3][4] * dz + S[1][4] * th1 + S[2][4] * th2 + S[0][4] * z) * (g * (l2 * l2) * (m2 * m2) * sin(th1) + J2 * g * m1 * sin(th1) * 4.0 + J2 * g * m2 * sin(th1) * 8.0 + g * (l2 * l2) * (m2 * m2) * sin(th1 - th2 * 2.0) - (dth2 * dth2) * (l2 * l2 * l2) * (m2 * m2) * sin(th1 - th2) + g * (l2 * l2) * m1 * m2 * sin(th1) - J2 * (dth2 * dth2) * l2 * m2 * sin(th1 - th2) * 4.0 + K[4] * dth1 * (l2 * l2) * (m2 * m2) * cos(th1) + K[5] * dth2 * (l2 * l2) * (m2 * m2) * cos(th1) + K[3] * dz * (l2 * l2) * (m2 * m2) * cos(th1) + K[1] * (l2 * l2) * (m2 * m2) * th1 * cos(th1) + K[2] * (l2 * l2) * (m2 * m2) * th2 * cos(th1) + K[0] * (l2 * l2) * (m2 * m2) * z * cos(th1) + J2 * K[4] * dth1 * m1 * cos(th1) * 4.0 + J2 * K[4] * dth1 * m2 * cos(th1) * 8.0 + J2 * K[5] * dth2 * m1 * cos(th1) * 4.0 + J2 * K[5] * dth2 * m2 * cos(th1) * 8.0 + J2 * K[3] * dz * m1 * cos(th1) * 4.0 + J2 * K[3] * dz * m2 * cos(th1) * 8.0 - (dth1 * dth1) * l1 * (l2 * l2) * (m2 * m2) * sin(th1 * 2.0 - th2 * 2.0) + J2 * K[1] * m1 * th1 * cos(th1) * 4.0 + J2 * K[1] * m2 * th1 * cos(th1) * 8.0 + J2 * K[2] * m1 * th2 * cos(th1) * 4.0 + J2 * K[2] * m2 * th2 * cos(th1) * 8.0 - K[4] * dth1 * (l2 * l2) * (m2 * m2) * cos(th1 - th2 * 2.0) - K[5] * dth2 * (l2 * l2) * (m2 * m2) * cos(th1 - th2 * 2.0) - K[3] * dz * (l2 * l2) * (m2 * m2) * cos(th1 - th2 * 2.0) + J2 * K[0] * m1 * z * cos(th1) * 4.0 + J2 * K[0] * m2 * z * cos(th1) * 8.0 - K[1] * (l2 * l2) * (m2 * m2) * th1 * cos(th1 - th2 * 2.0) - K[2] * (l2 * l2) * (m2 * m2) * th2 * cos(th1 - th2 * 2.0) - K[0] * (l2 * l2) * (m2 * m2) * z * cos(th1 - th2 * 2.0) + K[4] * dth1 * (l2 * l2) * m1 * m2 * cos(th1) + K[5] * dth2 * (l2 * l2) * m1 * m2 * cos(th1) + K[3] * dz * (l2 * l2) * m1 * m2 * cos(th1) + K[1] * (l2 * l2) * m1 * m2 * th1 * cos(th1) + K[2] * (l2 * l2) * m1 * m2 * th2 * cos(th1) + K[0] * (l2 * l2) * m1 * m2 * z * cos(th1)) * 2.0) / (J1 * J2 * 1.6E1 + (l1 * l1) * (l2 * l2) * (m2 * m2) * 2.0 + J2 * (l1 * l1) * m1 * 4.0 + J1 * (l2 * l2) * m2 * 4.0 + J2 * (l1 * l1) * m2 * 1.6E1 - (l1 * l1) * (l2 * l2) * (m2 * m2) * cos(th1 * 2.0 - th2 * 2.0) * 2.0 + (l1 * l1) * (l2 * l2) * m1 * m2) + (l2 * m2 * (S[4][5] * dth1 + S[5][5] * dth2 + S[3][5] * dz + S[1][5] * th1 + S[2][5] * th2 + S[0][5] * z) * (J1 * g * sin(th2) * 4.0 + J1 * K[4] * dth1 * cos(th2) * 4.0 + J1 * K[5] * dth2 * cos(th2) * 4.0 + J1 * K[3] * dz * cos(th2) * 4.0 + J1 * K[1] * th1 * cos(th2) * 4.0 + J1 * K[2] * th2 * cos(th2) * 4.0 + J1 * K[0] * z * cos(th2) * 4.0 + (dth1 * dth1) * (l1 * l1 * l1) * m1 * sin(th1 - th2) + (dth1 * dth1) * (l1 * l1 * l1) * m2 * sin(th1 - th2) * 4.0 - g * (l1 * l1) * m1 * sin(th1 * 2.0 - th2) - g * (l1 * l1) * m2 * sin(th1 * 2.0 - th2) * 2.0 + g * (l1 * l1) * m2 * sin(th2) * 2.0 + J1 * (dth1 * dth1) * l1 * sin(th1 - th2) * 4.0 + K[4] * dth1 * (l1 * l1) * m2 * cos(th2) * 2.0 + K[5] * dth2 * (l1 * l1) * m2 * cos(th2) * 2.0 + K[3] * dz * (l1 * l1) * m2 * cos(th2) * 2.0 + K[1] * (l1 * l1) * m2 * th1 * cos(th2) * 2.0 + K[2] * (l1 * l1) * m2 * th2 * cos(th2) * 2.0 + K[0] * (l1 * l1) * m2 * z * cos(th2) * 2.0 + (dth2 * dth2) * (l1 * l1) * l2 * m2 * sin(th1 * 2.0 - th2 * 2.0) - K[4] * dth1 * (l1 * l1) * m1 * cos(th1 * 2.0 - th2) - K[4] * dth1 * (l1 * l1) * m2 * cos(th1 * 2.0 - th2) * 2.0 - K[5] * dth2 * (l1 * l1) * m1 * cos(th1 * 2.0 - th2) - K[5] * dth2 * (l1 * l1) * m2 * cos(th1 * 2.0 - th2) * 2.0 - K[3] * dz * (l1 * l1) * m1 * cos(th1 * 2.0 - th2) - K[3] * dz * (l1 * l1) * m2 * cos(th1 * 2.0 - th2) * 2.0 - K[1] * (l1 * l1) * m1 * th1 * cos(th1 * 2.0 - th2) - K[1] * (l1 * l1) * m2 * th1 * cos(th1 * 2.0 - th2) * 2.0 - K[2] * (l1 * l1) * m1 * th2 * cos(th1 * 2.0 - th2) - K[2] * (l1 * l1) * m2 * th2 * cos(th1 * 2.0 - th2) * 2.0 - K[0] * (l1 * l1) * m1 * z * cos(th1 * 2.0 - th2) - K[0] * (l1 * l1) * m2 * z * cos(th1 * 2.0 - th2) * 2.0) * 2.0) / (J1 * J2 * 1.6E1 + (l1 * l1) * (l2 * l2) * (m2 * m2) * 2.0 + J2 * (l1 * l1) * m1 * 4.0 + J1 * (l2 * l2) * m2 * 4.0 + J2 * (l1 * l1) * m2 * 1.6E1 - (l1 * l1) * (l2 * l2) * (m2 * m2) * cos(th1 * 2.0 - th2 * 2.0) * 2.0 + (l1 * l1) * (l2 * l2) * m1 * m2);

    return dVxdt;
}

// some fundamental calculations
function scale_vec(s, v) {
    var sv = new Array(v.length);
    for (var k = 0; k < v.length; k++) {
        sv[k] = s * v[k];
    }
    return sv;
}

function add_vec(v1, v2) {
    var s = new Array(v1.length);
    for (var k = 0; k < v1.length; k++) {
        s[k] = v1[k] + v2[k];
    }
    return s;
}

function mult_vec(v1, v2) {
    var p = 0;
    for (var k = 0; k < v1.length; k++) {
        p = p + v1[k] * v2[k];
    }
    return p;
}

function quad_form(Q, x) {
    var s = 0;
    for (var i = 0; i < x.length; i++) {
        for (var j = 0; j < x.length; j++) {
            s = s + Q[i][j] * x[i] * x[j];
        }
    }
    return s;

}

function mod(n, m) {
    return ((n % m) + m) % m;
}

// one step 4th order Runge-Kutta routine 
function rk4(fc, x, t, dt) {
    let k1 = fc(t, x);
    let k2 = fc(t + dt / 2, add_vec(x, scale_vec(dt / 2, k1)));
    let k3 = fc(t + dt / 2, add_vec(x, scale_vec(dt / 2, k2)));
    let k4 = fc(t + dt, add_vec(x, scale_vec(dt, k3)));
    return add_vec(x, scale_vec(dt / 6, add_vec(add_vec(k1, scale_vec(2, k2)), add_vec(scale_vec(2, k3), k4))));
}

// controller
function wrap_x(x) {
    return [x[0] - r, mod(x[1] + PI, 2 * PI) - PI, mod(x[2] + PI, 2 * PI) - PI, x[3], x[4], x[5]];
}

function mode_transition(x, mode) {
    let new_mode = mode;
    let xs = wrap_x(x);
    Vt = 0.5 * quad_form(S, xs);
    dVtdt = dVdt(xs);

    Vt = 0.5 * quad_form(S, xs);
    dVtdt = dVdt(xs);
    switch (mode) {
        case 'lin':
            if (Vt > 20 && dVtdt > 0) {
                new_mode = 'swingup';
                mode_last_change = millis();
            }
            break;
        case 'swingup':
            if (Vt < 10 && dVtdt < 0 && cos(xs[1]) > 0.5 && cos(xs[2]) > 0.5) {
                new_mode = 'lin';
                mode_last_change = millis();
            }
            break;
    }
    return new_mode;
}

function u(r, x, mode) {
    Et = E(x);
    let xs = wrap_x(x);

    // controller
    switch (mode) {
        case 'lin':
            ut = -1 * mult_vec(K, xs);
            break;
        case 'swingup':
            let dExc = dEc(x);
            ut = -Et / dExc - 10 * xs[0] - 10 * xs[3];
            ut = constrain(ut, -10, 10);
            break;
    }
    // limit power
    ut = constrain(ut, -100, 100);
    return ut;
}

// visualization

function fmtnum(x) {
    return (x < 0 ? '' : ' ') + x.toFixed(3);
}

function draw_pend(x, alpha) {
    push();

    cw = 200;
    ch = 150;
    cr = 10;
    pw1 = 50;
    pw2 = 40;
    pr = 20;

    // draw cart
    translate(x[0] * 1000, -ch / 2);
    fill(232, 236, 239, 128 * alpha);
    stroke(77, 88, 89, 255 * alpha);
    rect(-cw / 2, -ch / 2, cw, ch, cr, cr, cr, cr);

    // draw pendulum
    rotate(x[1]);
    rect(-pw1 / 2, -pw1 / 2 - l1 * 1000, pw1, l1 * 1000 + pw1, pr, pr, pr, pr);
    translate(0, -l1 * 1000);
    rotate(x[2] - x[1]);
    rect(-pw2 / 2, -pw2 / 2 - l2 * 1000, pw2, l2 * 1000 + pw2, pr, pr, pr, pr);
    pop();
    push();
    // draw cart
    translate(x[0] * 1000, -ch / 2);
    noFill();
    stroke(77, 88, 89, 255 * alpha);
    strokeWeight(5);
    rect(-cw / 2, -ch / 2, cw, ch, cr, cr, cr, cr);
    fill(77, 88, 89, 128 * alpha);
    ellipse(0, 0, 30, 30);
    noFill();

    // draw pendulum
    rotate(x[1]);
    rect(-pw1 / 2, -pw1 / 2 - l1 * 1000, pw1, l1 * 1000 + pw1, pr, pr, pr, pr);
    translate(0, -l1 * 1000);
    rotate(x[2] - x[1]);
    rect(-pw2 / 2, -pw2 / 2 - l2 * 1000, pw2, l2 * 1000 + pw2, pr, pr, pr, pr);
    fill(77, 88, 89, 255 * alpha);
    ellipse(0, 0, 30, 30);
    noFill();
    pop();
}

function draw_note(x) {
    push();
    translate(x[0] * 1000, -ch / 2);
    textSize(48);
    fill(122, 126, 128);
    noStroke();
    textAlign(LEFT);
    let modestr = {
        'lin': 'Linear Quadratic Regulator',
        'swingup': 'Energy-based Control'
    }
    text(modestr[current_mode] + '\n' +
        'u(t)  = ' + fmtnum(ut / g) + ' g\n' +
        'E(t)  = ' + fmtnum(Et) + ' J\n' +
        'V(t)  = ' + Vt.toFixed(3) + '\n' +
        'dV(t) = ' + dVtdt.toFixed(3), cw / 2 + 50, -200);
    pop();
}

function set_view() {
    s = min(width / 3000, height / 3000);
    vymin = -height / s / 2;
    vymax = height / s / 2;
    xrng = width / s;
    yrng = height / s;
    let old_xrng = vxmax - vxmin;
    vxmin = vxmin - (xrng - old_xrng) / 2;
    vxmax = vxmax + (xrng - old_xrng) / 2;
}

function windowResized() {
    resizeCanvas(windowWidth, windowHeight);
    set_view();
}

var x_hist;
let hist_length = 20;

function setup() {
    var createRingBuffer = function (length) {

        var pointer = 0,
            buffer = [];

        return {
            get: function (key) {
                let idx = mod((pointer - key - 1), length);
                return buffer[idx];
            },
            push: function (item) {
                buffer[pointer] = item;
                pointer = (length + pointer + 1) % length;
            }
        };
    };
    x_hist = createRingBuffer(hist_length);
    for (var i = 0; i < hist_length; i++) {
        x_hist.push(xt.concat());
    }
    createCanvas(windowWidth, windowHeight);
    frameRate(60);
    set_view();
}

function draw_marker(x, t) {
    let ts = 50;
    textSize(54);
    if (x * 1e3 < vxmin) {
        triangle(vxmin, ts / 2, vxmin + ts / 2 * sqrt(3), 0, vxmin + ts / 2 * sqrt(3), ts);
        textAlign(LEFT);
        noStroke();
        text(t, vxmin, 48 + 40);
    } else if (x * 1e3 > vxmax) {
        triangle(vxmax, ts / 2, vxmax - ts / 2 * sqrt(3), 0, vxmax - ts / 2 * sqrt(3), ts);
        textAlign(RIGHT);
        noStroke();
        text(t, vxmax, 48 + 40);
    } else {
        triangle(x * 1e3, 0, x * 1e3 - ts / 2, +ts / 2 * sqrt(3), x * 1e3 + ts / 2, +ts / 2 * sqrt(3));
        textAlign(CENTER);
        noStroke();
        text(t, x * 1e3, 48 + 40);
    }
}

function draw_axis() {
    stroke(192, 198, 201);
    strokeWeight(2);
    strokeWeight(4);
    for (var gx = floor(vxmin / 1000) * 1000; gx <= floor(vxmax / 1000) * 1000; gx += 1000) {
        noStroke();
        stroke(122, 126, 128);
        line(gx, 0, gx, 40);
    }
    strokeWeight(4);
    line(vxmin, 0, vxmax, 0)
}

var r = 0;

function draw() {
    clear();
    resetMatrix();
    scale(s, s);

    // scroll viewport
    if (xt[0] * 1000 > vxmax - xrng * 0.25) {
        vxmax = xt[0] * 1000 + xrng * 0.25;
        vxmin = vxmax - xrng;
    }
    if (xt[0] * 1000 < vxmin + xrng * 0.25) {
        vxmin = xt[0] * 1000 - xrng * 0.25;
        vxmax = vxmin + xrng;
    }
    translate(-vxmin, -vymin);

    draw_axis();

    if (mouseIsPressed) {
        r = (mouseX / s + vxmin) / 1000;
        mode_last_change = millis();
    }

    if (current_mode === "lin" && mode_last_change + 5000 < millis()) {
        r = Math.random() * 20 - 10;
        mode_last_change = millis();
    }


    noStroke();
    fill(255, 0, 0, 128);
    draw_marker(r, '\nr(t) = ' + r.toFixed(3) + ' m');
    fill(0, 0, 255, 128);
    draw_marker(xt[0], 'x(t) = ' + xt[0].toFixed(3) + ' m');

    for (var i = hist_length - 1; i >= 0; i--) {
        draw_pend(x_hist.get(hist_length - i - 1), 0.1 / (hist_length + 1) * (i + 1) + 0.1);
    }
    draw_pend(xt, 1.0);
    draw_note(xt);
    x_hist.push(xt);
    current_mode = mode_transition(xt, current_mode);
    xt = rk4(function (t, x) {
        return add_vec(f(x), scale_vec(u(r, x, current_mode), B(x)));
    }, xt, millis() * 1e-3, 1 / 60);
}