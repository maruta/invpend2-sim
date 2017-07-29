let g = 9.8;
let m1 = 1;
let m2 = 1;
let l1 = 0.4;
let l2 = 0.6;
let J1 = 1 / 3 * m1 * l1 * l1 / 12;
let J2 = 1 / 3 * m2 * l2 * l2 / 12;


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

let K = [5.773503e+00, -1.245145e+02, 1.664657e+02, 7.810893e+00, -2.535013e+00, 2.379305e+01, ];
let S = [
    [1.352886e+00, -4.390771e-01, 4.121077e+00, 9.146507e-01, 3.943857e-01, 7.367118e-01, ],
    [-4.390771e-01, 1.716117e+01, -2.483344e+01, -9.884071e-01, 4.639436e-02, -3.637064e+00, ],
    [4.121077e+00, -2.483344e+01, 5.356611e+01, 4.838637e+00, 1.827593e+00, 8.441639e+00, ],
    [9.146507e-01, -9.884071e-01, 4.838637e+00, 8.617041e-01, 3.812820e-01, 8.374408e-01, ],
    [3.943857e-01, 4.639436e-02, 1.827593e+00, 3.812820e-01, 2.178164e-01, 3.399085e-01, ],
    [7.367118e-01, -3.637064e+00, 8.441639e+00, 8.374408e-01, 3.399085e-01, 1.356479e+00, ],
];

function dVdt(x) {
    let z = x[0];
    let th1 = x[1];
    let th2 = x[2];
    let dz = x[3];
    let dth1 = x[4];
    let dth2 = x[5];

    let dVxdt = dth2 * (dth1 * 1.827593419238315 + dth2 * 8.441638806461018 + dz * 4.838637263239541 - th1 * 2.483344021304082E1 + th2 * 5.356611299640473E1 + z * 4.121077239378319) + dz * (dth1 * 3.943856994488532E-1 + dth2 * 7.367117758101513E-1 + dz * 9.146507034066623E-1 - th1 * 4.390771154279364E-1 + th2 * 4.121077239378319 + z * 1.352886324423944) - (dth1 * (-2.535012907873332) + dth2 * 2.379305053506306E1 + dz * 7.810892835891162 - th1 * 1.245144929728302E2 + th2 * 1.664657282958085E2 + z * 5.773502691896171) * (dth1 * 3.812820205527619E-1 + dth2 * 8.374408097330027E-1 + dz * 8.617041299013249E-1 - th1 * 9.884071242788216E-1 + th2 * 4.838637263239541 + z * 9.146507034066623E-1) - dth1 * (dth1 * -4.639435664626901E-2 + dth2 * 3.637064125960309 + dz * 9.884071242788216E-1 - th1 * 1.716117492535966E1 + th2 * 2.483344021304082E1 + z * 4.390771154279364E-1) + (l1 * (dth1 * 2.178164076428267E-1 + dth2 * 3.399084634881735E-1 + dz * 3.812820205527619E-1 + th1 * 4.639435664626901E-2 + th2 * 1.827593419238315 + z * 3.943856994488532E-1) * (dth1 * (l2 * l2) * (m2 * m2) * cos(th1) * -5.708341593638867E15 + dth2 * (l2 * l2) * (m2 * m2) * cos(th1) * 5.35771867618587E16 + dz * (l2 * l2) * (m2 * m2) * cos(th1) * 1.758856703257516E16 - (l2 * l2) * (m2 * m2) * th1 * cos(th1) * 2.803817120773322E17 + (l2 * l2) * (m2 * m2) * th2 * cos(th1) * 3.748474959614806E17 + g * (l2 * l2) * (m2 * m2) * sin(th1) * 2.251799813685248E15 + (l2 * l2) * (m2 * m2) * z * cos(th1) * 1.300077228592308E16 - J2 * dth1 * m1 * cos(th1) * 2.283336637455547E16 - J2 * dth1 * m2 * cos(th1) * 4.566673274911094E16 + J2 * dth2 * m1 * cos(th1) * 2.143087470474348E17 + J2 * dth2 * m2 * cos(th1) * 4.286174940948696E17 + J2 * dz * m1 * cos(th1) * 7.035426813030063E16 + J2 * dz * m2 * cos(th1) * 1.407085362606013E17 - J2 * m1 * th1 * cos(th1) * 1.121526848309329E18 + J2 * m1 * th2 * cos(th1) * 1.499389983845923E18 - J2 * m2 * th1 * cos(th1) * 2.243053696618658E18 + J2 * m2 * th2 * cos(th1) * 2.998779967691845E18 + J2 * g * m1 * sin(th1) * 9.007199254740992E15 + J2 * g * m2 * sin(th1) * 1.801439850948198E16 + dth1 * (l2 * l2) * (m2 * m2) * cos(th1 - th2 * 2.0) * 5.708341593638867E15 - dth2 * (l2 * l2) * (m2 * m2) * cos(th1 - th2 * 2.0) * 5.35771867618587E16 - dz * (l2 * l2) * (m2 * m2) * cos(th1 - th2 * 2.0) * 1.758856703257516E16 + J2 * m1 * z * cos(th1) * 5.20030891436923E16 + J2 * m2 * z * cos(th1) * 1.040061782873846E17 + (l2 * l2) * (m2 * m2) * th1 * cos(th1 - th2 * 2.0) * 2.803817120773322E17 - (l2 * l2) * (m2 * m2) * th2 * cos(th1 - th2 * 2.0) * 3.748474959614806E17 + g * (l2 * l2) * (m2 * m2) * sin(th1 - th2 * 2.0) * 2.251799813685248E15 - (l2 * l2) * (m2 * m2) * z * cos(th1 - th2 * 2.0) * 1.300077228592308E16 - (dth2 * dth2) * (l2 * l2 * l2) * (m2 * m2) * sin(th1 - th2) * 2.251799813685248E15 - dth1 * (l2 * l2) * m1 * m2 * cos(th1) * 5.708341593638867E15 + dth2 * (l2 * l2) * m1 * m2 * cos(th1) * 5.35771867618587E16 + dz * (l2 * l2) * m1 * m2 * cos(th1) * 1.758856703257516E16 - (l2 * l2) * m1 * m2 * th1 * cos(th1) * 2.803817120773322E17 + (l2 * l2) * m1 * m2 * th2 * cos(th1) * 3.748474959614806E17 + g * (l2 * l2) * m1 * m2 * sin(th1) * 2.251799813685248E15 + (l2 * l2) * m1 * m2 * z * cos(th1) * 1.300077228592308E16 - J2 * (dth2 * dth2) * l2 * m2 * sin(th1 - th2) * 9.007199254740992E15 - (dth1 * dth1) * l1 * (l2 * l2) * (m2 * m2) * sin(th1 * 2.0 - th2 * 2.0) * 2.251799813685248E15) * 8.881784197001252E-16) / (J1 * J2 * 1.6E1 + (l1 * l1) * (l2 * l2) * (m2 * m2) * 2.0 + J2 * (l1 * l1) * m1 * 4.0 + J1 * (l2 * l2) * m2 * 4.0 + J2 * (l1 * l1) * m2 * 1.6E1 - (l1 * l1) * (l2 * l2) * (m2 * m2) * cos(th1 * 2.0 - th2 * 2.0) * 2.0 + (l1 * l1) * (l2 * l2) * m1 * m2) + (l2 * m2 * (dth1 * 3.399084634881735E-1 + dth2 * 1.356479074543226 + dz * 8.374408097330027E-1 - th1 * 3.637064125960309 + th2 * 8.441638806461018 + z * 7.367117758101513E-1) * (J1 * dth1 * cos(th2) * -2.283336637455547E16 + J1 * dth2 * cos(th2) * 2.143087470474348E17 + J1 * dz * cos(th2) * 7.035426813030063E16 - J1 * th1 * cos(th2) * 1.121526848309329E18 + J1 * th2 * cos(th2) * 1.499389983845923E18 + J1 * g * sin(th2) * 9.007199254740992E15 + J1 * z * cos(th2) * 5.20030891436923E16 + dth1 * (l1 * l1) * m1 * cos(th1 * 2.0 - th2) * 5.708341593638867E15 + dth1 * (l1 * l1) * m2 * cos(th1 * 2.0 - th2) * 1.141668318727773E16 - dth2 * (l1 * l1) * m1 * cos(th1 * 2.0 - th2) * 5.35771867618587E16 - dth2 * (l1 * l1) * m2 * cos(th1 * 2.0 - th2) * 1.071543735237174E17 - dz * (l1 * l1) * m1 * cos(th1 * 2.0 - th2) * 1.758856703257516E16 - dz * (l1 * l1) * m2 * cos(th1 * 2.0 - th2) * 3.517713406515032E16 + (l1 * l1) * m1 * th1 * cos(th1 * 2.0 - th2) * 2.803817120773322E17 - (l1 * l1) * m1 * th2 * cos(th1 * 2.0 - th2) * 3.748474959614806E17 + (l1 * l1) * m2 * th1 * cos(th1 * 2.0 - th2) * 5.607634241546644E17 - (l1 * l1) * m2 * th2 * cos(th1 * 2.0 - th2) * 7.496949919229613E17 + (dth1 * dth1) * (l1 * l1 * l1) * m1 * sin(th1 - th2) * 2.251799813685248E15 + (dth1 * dth1) * (l1 * l1 * l1) * m2 * sin(th1 - th2) * 9.007199254740992E15 - g * (l1 * l1) * m1 * sin(th1 * 2.0 - th2) * 2.251799813685248E15 - g * (l1 * l1) * m2 * sin(th1 * 2.0 - th2) * 4.503599627370496E15 - (l1 * l1) * m1 * z * cos(th1 * 2.0 - th2) * 1.300077228592308E16 - (l1 * l1) * m2 * z * cos(th1 * 2.0 - th2) * 2.600154457184615E16 - dth1 * (l1 * l1) * m2 * cos(th2) * 1.141668318727773E16 + dth2 * (l1 * l1) * m2 * cos(th2) * 1.071543735237174E17 + dz * (l1 * l1) * m2 * cos(th2) * 3.517713406515032E16 - (l1 * l1) * m2 * th1 * cos(th2) * 5.607634241546644E17 + (l1 * l1) * m2 * th2 * cos(th2) * 7.496949919229613E17 + g * (l1 * l1) * m2 * sin(th2) * 4.503599627370496E15 + (l1 * l1) * m2 * z * cos(th2) * 2.600154457184615E16 + J1 * (dth1 * dth1) * l1 * sin(th1 - th2) * 9.007199254740992E15 + (dth2 * dth2) * (l1 * l1) * l2 * m2 * sin(th1 * 2.0 - th2 * 2.0) * 2.251799813685248E15) * 8.881784197001252E-16) / (J1 * J2 * 1.6E1 + (l1 * l1) * (l2 * l2) * (m2 * m2) * 2.0 + J2 * (l1 * l1) * m1 * 4.0 + J1 * (l2 * l2) * m2 * 4.0 + J2 * (l1 * l1) * m2 * 1.6E1 - (l1 * l1) * (l2 * l2) * (m2 * m2) * cos(th1 * 2.0 - th2 * 2.0) * 2.0 + (l1 * l1) * (l2 * l2) * m1 * m2);

    return dVxdt;
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
var mode = 'swingup';
var mode_last_change = 0;

function mod(n, m) {
    return ((n % m) + m) % m;
}

var Ex = 0;
var Vx = 0;
var ux = 0;
var dVxdt = 0;

function u(r, x) {
    Ex = E(x);
    let dExc = dEc(x);
    let xs = [x[0] - r, mod(x[1] + PI, 2 * PI) - PI, mod(x[2] + PI, 2 * PI) - PI, x[3], x[4], x[5]];

    Vx = 0.5 * quad_form(S, xs);
    dVxdt = dVdt(xs);
    // mode transition rule
    switch (mode) {
        case 'lin':
            if (dVxdt > 1e-5 || Vx > 10) {
                mode = 'swingup';
                mode_last_change = millis();
            }
            break;
        case 'swingup':
            if (Vx < 5 && dVxdt < -1e-5 && cos(xs[1]) > 0 && cos(xs[2]) > 0) {
                mode = 'lin';
                mode_last_change = millis();
            }
            break;
    }
    // controller
    switch (mode) {
        case 'lin':
            ux = -1 * mult_vec(K, xs);
            break;
        case 'swingup':
            ux = -Ex / dExc - 10 * xs[0] - 10 * xs[3];
            ux = constrain(ux, -10, 10);
            break;
    }
    // limit power
    ux = constrain(ux, -1000, 1000);
    return ux;
}

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

function rk4(fc, x, t, dt) {
    let k1 = fc(t, x);
    let k2 = fc(t + dt / 2, add_vec(x, scale_vec(dt / 2, k1)));
    let k3 = fc(t + dt / 2, add_vec(x, scale_vec(dt / 2, k2)));
    let k4 = fc(t + dt, add_vec(x, scale_vec(dt, k3)));
    return add_vec(x, scale_vec(dt / 6, add_vec(add_vec(k1, scale_vec(2, k2)), add_vec(scale_vec(2, k3), k4))));
}

function fmtnum(x) {
    return (x < 0 ? '' : ' ') + x.toFixed(3);
}

function draw_pend(x) {
    push();

    cw = 200;
    ch = 150;
    cr = 10;
    pw1 = 50;
    pw2 = 40;
    pr = 20;

    // draw cart
    translate(x[0] * 1000, -ch / 2);
    fill(232, 236, 239, 128);
    stroke(77, 88, 89);
    rect(-cw / 2, -ch / 2, cw, ch, cr, cr, cr, cr);

    // draw pendulum
    rotate(x[1]);
    rect(-pw1 / 2, -pw1 / 2 - l1 * 1000, pw1, l1 * 1000 + pw1, pr, pr, pr, pr);
    translate(0, -l1 * 1000);
    rotate(x[2] - x[1]);
    rect(-pw2 / 2, -pw2 / 2 - l2 * 1000, pw2, l2 * 1000 + pw2, pr, pr, pr, pr);
    pop();

    // draw cart
    translate(x[0] * 1000, -ch / 2);
    noFill();
    stroke(77, 88, 89);
    strokeWeight(5);
    rect(-cw / 2, -ch / 2, cw, ch, cr, cr, cr, cr);
    fill(77, 88, 89, 128);
    ellipse(0, 0, 30, 30);
    noFill();
    push();
    textSize(48);
    fill(122, 126, 128);
    noStroke();
    textAlign(LEFT);
    let modestr = {
        'lin': 'Linear Quadratic Regulator',
        'swingup': 'Energy-based Control'
    }
    text(modestr[mode] + '\n' +
        'u(t)  = ' + fmtnum(ux / g) + ' g\n' +
        'E(t)  = ' + fmtnum(Ex) + ' J\n' +
        'V(t)  = ' + fmtnum(Vx) + '\n' +
        'dV(t) = ' + fmtnum(dVxdt), cw / 2 + 50, -200);
    pop();
    // draw pendulum
    rotate(x[1]);
    rect(-pw1 / 2, -pw1 / 2 - l1 * 1000, pw1, l1 * 1000 + pw1, pr, pr, pr, pr);
    translate(0, -l1 * 1000);
    rotate(x[2] - x[1]);
    rect(-pw2 / 2, -pw2 / 2 - l2 * 1000, pw2, l2 * 1000 + pw2, pr, pr, pr, pr);
    fill(77, 88, 89);
    ellipse(0, 0, 30, 30);
    noFill();
    pop();


}

var x = [0, 2 * 3.14 * Math.random(), 2 * 3.14 * Math.random(), 0, 0, 0];
var xmin, xmax, ymin, ymax;

function preload() {}

function windowResized() {
    resizeCanvas(windowWidth, windowHeight);
    s = min(width / 3000, height / 3000);
    ymin = -height / s / 2;
    ymax = height / s / 2;
}

function setup() {
    //    textFont('Share Tech Mono');
    createCanvas(windowWidth, windowHeight);
    frameRate(60);
    s = min(width / 3000, height / 3000);
    xmin = -width / s / 2;
    xmax = width / s / 2;
    ymin = -height / s / 2;
    ymax = height / s / 2;
}


var r = 0;

function draw() {
    clear();
    resetMatrix();
    s = min(width / 3000, height / 3000);
    scale(s, s);

    xrng = width / s;
    yrng = height / s;
    if (x[0] * 1000 > xmax - xrng * 0.25) {
        xmax = x[0] * 1000 + xrng * 0.25;
        xmin = xmax - xrng;
    }
    if (x[0] * 1000 < xmin + xrng * 0.25) {
        xmin = x[0] * 1000 - xrng * 0.25;
        xmax = xmin + xrng;
    }
    translate(-xmin, -ymin);
    stroke(192, 198, 201);
    strokeWeight(2);
    textSize(54);
    fill(122, 126, 128);
    strokeWeight(5);
    textAlign(CENTER);
    for (var gx = floor(xmin / 1000) * 1000; gx <= floor(xmax / 1000) * 1000; gx += 1000) {
        noStroke();
        text(gx / 1000 + ' m', gx, +48 + 40);
        stroke(122, 126, 128);
        line(gx, 0, gx, 20);
    }
    strokeWeight(3);
    line(xmin, 0, xmax, 0)

    if (mouseIsPressed) {
        r = (mouseX / s + xmin) / 1000;
        mode_last_change = millis();
    }

    if (mode === "lin" && mode_last_change + 5000 < millis()) {
        r = Math.random() * 20 - 10;
        mode_last_change = millis();
    }

    let ts = 50;
    noStroke();
    fill(255, 0, 0, 128);
    triangle(r * 1e3, 0, r * 1e3 - ts / 2, +ts / 2 * sqrt(3), r * 1e3 + ts / 2, +ts / 2 * sqrt(3));
    fill(0, 0, 255, 128);
    triangle(x[0] * 1e3, 0, x[0] * 1e3 - ts / 2, +ts / 2 * sqrt(3), x[0] * 1e3 + ts / 2, +ts / 2 * sqrt(3));
    draw_pend(x);
    x = rk4(function (t, x) {
        return add_vec(f(x), scale_vec(u(r, x), B(x)));
    }, x, 0, 1 / 60);
}