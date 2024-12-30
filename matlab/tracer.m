RES = [32 32];
FOVY = deg2rad(90);
EYE = [-4;2;4];
TARGET = [0;0;0];
CLIP = [0.1 10.050000000745058];

function val = torus(p)
    % torus
    R = 2;
    r = 0.5;
    x = p(1); y = p(2); z = p(3);
    lhs = 4*R*R*(x*x+y*y);
    rhs = (x*x+y*y+z*z+R*R-r*r);
    val = rhs*rhs-lhs;
end


function val = barth_sextic(p)
    % torus
    x = p(1);
    y = p(2);
    z = p(3);
    val = 4 * (2.618 * x * x - y * x) * (2.618 * y * y - z * z) * (2.618 * z * z - x * x) - 4.236 * power(x * x + y * y + z * z - 1, 2);
end

function val = sphere(p)
    R = 1;
    val = norm(p)-R;
end


function val = surface_func(p)
    val = barth_sextic(p);
end



f = cot(FOVY/2);
aspect = RES(1)/RES(2);
near = CLIP(1);
far = CLIP(2);
perspective = [f/aspect 0 0 0; 0 f 0 0; 0 0 (far+near)/(near-far) (2*far*near)/((near-far)); 0 0 -1 0];

fwd = EYE-TARGET;
fwd = fwd / norm(fwd);
right = cross(fwd,[0;1;0]);
right = right/norm(right);
up = cross(right, fwd);
up = up/norm(up);
rot = [right,up,fwd];
view = zeros(4);
view(1:3,1:3) = transpose(rot);
view(1:3,4) = -EYE;
view(4,4) = 1;


yedge = tan(FOVY / 2);
xedge = RES(1)/RES(2)*yedge;
invView = transpose([0.707107, -0, 0.707107, -0; 0.235702, 0.942809, -0.235702, 0; -0.666667, 0.333333, 0.666667, -0; -4, 2, 4, 1]);

viewProj = perspective * view;
invViewProj = inv(viewProj);

img = zeros([RES(1),RES(2),3]);

for Y = 1:RES(2)
    for X = 1:RES(1)
        %if Y ~= 15 || X ~= 15
        %    continue
        %end

       ndc = zeros(4,1);
       ndc(1) = (X-1)/(RES(1)-1)*2-1;
       ndc(2) = 1-2*(Y-1)/(RES(2)-1);
       ndc(4) = 1;

       wp = invViewProj * ndc;
       wp = wp/wp(4);
       wp = wp(1:3);

       ahead = [ndc(1)*xedge; ndc(2)*yedge; -1; 1];
       wp = invView * ahead;
       wp = wp(1:3,1);

       rayDir = wp-EYE;
       rayDir = rayDir/norm(rayDir);

       rayFunc = @(t) surface_func(EYE+rayDir*t);
       %chebfunpref.setDefaults({'techPrefs','eps'}, 1e-15);
       approx = chebfun(rayFunc, [near, far], 'splitting', 'on');


       seg_count = size(approx.funs);
       seg_count = seg_count(1);
       for i=1:seg_count
           f = approx.funs(i,1);
           degree = size(f{1,1}.onefun.coeffs); %1 or 2?
           degree = degree(1)-1;
           next_cheb_points = cos(linspace(0,1,degree)*pi);
           next_cheb_points = cos(linspace(0,1,degree*2-1)*pi);

           next_cheb_points = (next_cheb_points/2+0.5)*(far-near)+near;
           rayFunc(next_cheb_points);
           approx(next_cheb_points);

           errors = zeros(size(next_cheb_points));
           for e = 1:size(errors,2)
            errors(e) = abs(rayFunc(next_cheb_points(e))-approx(next_cheb_points(e)));
           end
           err = norm(errors, inf)
       end

       r = roots(approx)
       if size(r,1) > 0
        dist = r(1);
        wp = EYE+rayDir*dist;
        n = zeros(3,1);
        n(1) = surface_func([wp(1) + 0.01; wp(2); wp(3)]) - surface_func([wp(1) - 0.01; wp(2); wp(3)]);
        n(2) = surface_func([wp(1); wp(2)+ 0.01; wp(3)]) - surface_func([wp(1); wp(2) - 0.01; wp(3)]);
        n(3) = surface_func([wp(1); wp(2); wp(3) + 0.01]) - surface_func([wp(1); wp(2); wp(3) - 0.01]);
        n = n/norm(n,2);

        to_eye = EYE-wp;
        to_eye = to_eye/norm(to_eye,2);
        if dot(n, to_eye) < 0
            n = -1*n;
        end

          LIGHT_DIR = [0.08105321228504181; -0.7479685544967651;-0.6379234790802002];
          LIGHT_DIR = LIGHT_DIR/norm(LIGHT_DIR,2);
          diffuse = dot(n, -LIGHT_DIR);
            
          img(Y,X,1) = diffuse;
          img(Y,X,2) = diffuse;
          img(Y,X,3) = diffuse;
        
       else
          img(Y,X,2) = 1;
       end
    end
end

imwrite(img,"asd_cheb.png")