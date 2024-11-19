for i = 2:nx-1
    for j = 2:ny-1
        for k = 2:nz-1
            ux(i,j,k) = sum(f(i,j,k,[2,16,10,8,14])) - sum(f(i,j,k,[3,11,17,15,8]));
            uy(i,j,k) = sum(f(i,j,k,[4,8,15,18,12])) - sum(f(i,j,k,[5,14,9,13,19]));
            uz(i,j,k) = sum(f(i,j,k,[7,16,11,18,13])) - sum(f(i,j,k,[6,10,17,12,19]));
            ux(i,j,k) = ux(i,j,k) ./ rho(i,j,k) + ffx(i,j,k) * 0.5 ./ rho(i,j,k);
            uy(i,j,k) = uy(i,j,k) ./ rho(i,j,k) + ffy(i,j,k) * 0.5 ./ rho(i,j,k);
            uz(i,j,k) = uz(i,j,k) ./ rho(i,j,k) + ffz(i,j,k) * 0.5 ./ rho(i,j,k);
            uu = 0.5 * (ux(i,j,k).^2 + uy(i,j,k).^2 + uz(i,j,k).^2) / cssq;
            rho(i,j,k) = sum(f(i,j,k,:),4);
            for l = 1:fpoints
                udotc = (ux(i,j,k) * cix(l) + uy(i,j,k) * ciy(l) + uz(i,j,k) * ciz(l)) / cssq;
                HeF = 0.5 .* (w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu))) ...
                    .* ((cix(l) - ux(i,j,k)) .* ffx(i,j,k) + ...
                        (ciy(l) - uy(i,j,k)) .* ffy(i,j,k) + ...
                        (ciz(l) - uz(i,j,k)) .* ffz(i,j,k) ...
                       ) ./ (rho(i,j,k) .* cssq);
                feq = w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu)) - HeF;
                fneq(l) = f(i,j,k,l) - feq;
            end
            pxx(i,j,k) = sum(fneq([2,3,8,9,10,11,14,15,16,17]));
            pyy(i,j,k) = sum(fneq([4,5,8,9,12,13,14,15,18,19]));
            pzz(i,j,k) = sum(fneq([6,7,10,11,12,13,16,17,18,19]));
            pxy(i,j,k) = fneq(8) + fneq(9) - fneq(14) - fneq(15);
            pxz(i,j,k) = fneq(10) + fneq(11) - fneq(16) - fneq(17);
            pyz(i,j,k) = fneq(12) + fneq(13) - fneq(18) - fneq(19);
        end
    end
end

ux(2:nx-1,2:ny-1,2:nz-1) = sum(f(2:nx-1,2:ny-1,2:nz-1,[2,16,10,8,14])) - sum(f(2:nx-1,2:ny-1,2:nz-1,[3,11,17,15,8]));
uy(2:nx-1,2:ny-1,2:nz-1) = sum(f(2:nx-1,2:ny-1,2:nz-1,[4,8,15,18,12])) - sum(f(2:nx-1,2:ny-1,2:nz-1,[5,14,9,13,19]));
uz(2:nx-1,2:ny-1,2:nz-1) = sum(f(2:nx-1,2:ny-1,2:nz-1,[7,16,11,18,13])) - sum(f(2:nx-1,2:ny-1,2:nz-1,[6,10,17,12,19]));
ux(2:nx-1,2:ny-1,2:nz-1) = ux(2:nx-1,2:ny-1,2:nz-1) ./ rho(2:nx-1,2:ny-1,2:nz-1) + ffx(2:nx-1,2:ny-1,2:nz-1) * 0.5 ./ rho(2:nx-1,2:ny-1,2:nz-1);
uy(2:nx-1,2:ny-1,2:nz-1) = uy(2:nx-1,2:ny-1,2:nz-1) ./ rho(2:nx-1,2:ny-1,2:nz-1) + ffy(2:nx-1,2:ny-1,2:nz-1) * 0.5 ./ rho(2:nx-1,2:ny-1,2:nz-1);
uz(2:nx-1,2:ny-1,2:nz-1) = uz(2:nx-1,2:ny-1,2:nz-1) ./ rho(2:nx-1,2:ny-1,2:nz-1) + ffz(2:nx-1,2:ny-1,2:nz-1) * 0.5 ./ rho(2:nx-1,2:ny-1,2:nz-1);
uu = 0.5 * (ux(2:nx-1,2:ny-1,2:nz-1).^2 + uy(2:nx-1,2:ny-1,2:nz-1).^2 + uz(2:nx-1,2:ny-1,2:nz-1).^2) / cssq;
rho(2:nx-1,2:ny-1,2:nz-1) = sum(f(2:nx-1,2:ny-1,2:nz-1,:),4);
for l = 1:fpoints
    udotc = (ux(2:nx-1,2:ny-1,2:nz-1) * cix(l) + uy(2:nx-1,2:ny-1,2:nz-1) * ciy(l) + uz(2:nx-1,2:ny-1,2:nz-1) * ciz(l)) / cssq;
    HeF = 0.5 .* (w(l) * (rho(2:nx-1,2:ny-1,2:nz-1) + rho(2:nx-1,2:ny-1,2:nz-1) .* (udotc + 0.5.*udotc.^2 - uu))) ...
        .* ((cix(l) - ux(2:nx-1,2:ny-1,2:nz-1)) .* ffx(2:nx-1,2:ny-1,2:nz-1) + ...
            (ciy(l) - uy(2:nx-1,2:ny-1,2:nz-1)) .* ffy(2:nx-1,2:ny-1,2:nz-1) + ...
            (ciz(l) - uz(2:nx-1,2:ny-1,2:nz-1)) .* ffz(2:nx-1,2:ny-1,2:nz-1) ...
           ) ./ (rho(2:nx-1,2:ny-1,2:nz-1) .* cssq);
    feq = w(l) * (rho(2:nx-1,2:ny-1,2:nz-1) + rho(2:nx-1,2:ny-1,2:nz-1) .* (udotc + 0.5.*udotc.^2 - uu)) - HeF;
    fneq(l) = f(2:nx-1,2:ny-1,2:nz-1,l) - feq;
end
pxx(2:nx-1,2:ny-1,2:nz-1) = sum(fneq([2,3,8,9,10,11,14,15,16,17]));
pyy(2:nx-1,2:ny-1,2:nz-1) = sum(fneq([4,5,8,9,12,13,14,15,18,19]));
pzz(2:nx-1,2:ny-1,2:nz-1) = sum(fneq([6,7,10,11,12,13,16,17,18,19]));
pxy(2:nx-1,2:ny-1,2:nz-1) = fneq(8) + fneq(9) - fneq(14) - fneq(15);
pxz(2:nx-1,2:ny-1,2:nz-1) = fneq(10) + fneq(11) - fneq(16) - fneq(17);
pyz(2:nx-1,2:ny-1,2:nz-1) = fneq(12) + fneq(13) - fneq(18) - fneq(19);

%% 

for i = 2:nx-1
    for j = 2:ny-1
        for k = 2:nz-1
            uu = 0.5 * (ux(i,j,k).^2 + uy(i,j,k).^2 + uz(i,j,k).^2) / cssq;
            for l = 1:fpoints
                udotc = (ux(i,j,k) * cix(l) + uy(i,j,k) * ciy(l) + uz(i,j,k) * ciz(l)) / cssq;
                feq = w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu));
                HeF = 0.5 .* (w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu))) .* ...
                    ((cix(l) - ux(i,j,k)) .* ffx(i,j,k) + ...
                     (ciy(l) - uy(i,j,k)) .* ffy(i,j,k) + ...
                     (ciz(l) - uz(i,j,k)) .* ffz(i,j,k) ...
                    ) ./ (rho(i,j,k) .* cssq);
                fneq = (cix(l) .* cix(l) - cssq) * pxx(i,j,k) + ...
                       (ciy(l) .* ciy(l) - cssq) * pyy(i,j,k) + ...
                       (ciz(l) .* ciz(l) - cssq) * pzz(i,j,k) + ...
                       (cix(l) .* ciy(l) - cssq) * pxy(i,j,k) + ...
                       (cix(l) .* ciz(l) - cssq) * pxz(i,j,k) + ...
                       (ciy(l) .* ciz(l) - cssq) * pyz(i,j,k);
                f(i+cix(l),j+ciy(l),k+ciz(l),l) = feq + (1-omega) * (w(l) / (2*cssq^2)) * fneq + HeF;
            end
            for l = 1:gpoints
                udotc = (ux(i,j,k) * cix(l) + uy(i,j,k) * ciy(l) + uz(i,j,k) * ciz(l)) / cssq;
                feq = w_g(l) .* phi(i,j,k) .* (1 + udotc);
                Hi = sharp_c .* phi(i,j,k) .* (1 - phi(i,j,k)) .* (cix(l) .* normx(i,j,k) + ciy(l) .* normy(i,j,k) + ciz(l) .* normz(i,j,k)); 
                g(i,j,k,l) = feq + w_g(l) .* Hi;
            end
        end
    end
end

uu = 0.5 * (ux(2:nx-1,2:ny-1,2:nz-1).^2 + uy(2:nx-1,2:ny-1,2:nz-1).^2 + uz(2:nx-1,2:ny-1,2:nz-1).^2) / cssq;
for l = 1:fpoints
    udotc = (ux(2:nx-1,2:ny-1,2:nz-1) * cix(l) + uy(2:nx-1,2:ny-1,2:nz-1) * ciy(l) + uz(2:nx-1,2:ny-1,2:nz-1) * ciz(l)) / cssq;
    feq = w(l) * (rho(2:nx-1,2:ny-1,2:nz-1) + rho(2:nx-1,2:ny-1,2:nz-1) .* (udotc + 0.5.*udotc.^2 - uu));
    HeF = 0.5 .* (w(l) * (rho(2:nx-1,2:ny-1,2:nz-1) + rho(2:nx-1,2:ny-1,2:nz-1) .* (udotc + 0.5.*udotc.^2 - uu))) .* ...
        ((cix(l) - ux(2:nx-1,2:ny-1,2:nz-1)) .* ffx(2:nx-1,2:ny-1,2:nz-1) + ...
         (ciy(l) - uy(2:nx-1,2:ny-1,2:nz-1)) .* ffy(2:nx-1,2:ny-1,2:nz-1) + ...
         (ciz(l) - uz(2:nx-1,2:ny-1,2:nz-1)) .* ffz(2:nx-1,2:ny-1,2:nz-1) ...
        ) ./ (rho(2:nx-1,2:ny-1,2:nz-1) .* cssq);
    fneq = (cix(l) .* cix(l) - cssq) * pxx(2:nx-1,2:ny-1,2:nz-1) + ...
           (ciy(l) .* ciy(l) - cssq) * pyy(2:nx-1,2:ny-1,2:nz-1) + ...
           (ciz(l) .* ciz(l) - cssq) * pzz(2:nx-1,2:ny-1,2:nz-1) + ...
           (cix(l) .* ciy(l) - cssq) * pxy(2:nx-1,2:ny-1,2:nz-1) + ...
           (cix(l) .* ciz(l) - cssq) * pxz(2:nx-1,2:ny-1,2:nz-1) + ...
           (ciy(l) .* ciz(l) - cssq) * pyz(2:nx-1,2:ny-1,2:nz-1);
    f((2:nx-1)+cix(l),(2:ny-1)+ciy(l),(2:nz-1)+ciz(l),l) = feq + (1-omega) * (w(l) / (2*cssq^2)) * fneq + HeF;
end
for l = 1:gpoints
    udotc = (ux(2:nx-1,2:ny-1,2:nz-1) * cix(l) + uy(2:nx-1,2:ny-1,2:nz-1) * ciy(l) + uz(2:nx-1,2:ny-1,2:nz-1) * ciz(l)) / cssq;
    feq = w_g(l) .* phi(2:nx-1,2:ny-1,2:nz-1) .* (1 + udotc);
    Hi = sharp_c .* phi(2:nx-1,2:ny-1,2:nz-1) .* (1 - phi(2:nx-1,2:ny-1,2:nz-1)) .* (cix(l) .* normx(2:nx-1,2:ny-1,2:nz-1) + ciy(l) .* normy(2:nx-1,2:ny-1,2:nz-1) + ciz(l) .* normz(2:nx-1,2:ny-1,2:nz-1)); 
    g(2:nx-1,2:ny-1,2:nz-1,l) = feq + w_g(l) .* Hi;
end