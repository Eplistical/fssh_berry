% exact solution for 2D model

function exact_2d_tully1(r, k, sigma, init_s, mass, fID)
    tic;

    xI = r(1);
    yI = r(2);
    kxI = k(1);
    kyI = k(2);
    sigmax = sigma(1);
    sigmay = sigma(2);
    c1 = sqrt(1 - init_s);
    c2 = sqrt(init_s);

    N = 2;
    L = 16;
    M = 256;                    

    dt = 0.1;
    Nstep = ceil(7.5 / (k(1) / mass) / dt); % 5000;
    tgraph = 100;

    enable_pop_only = false;
    enable_plot = true;
    enable_plot_diab_surf = false;
    enable_plot_surf = false;
    % grids
    x0 = linspace(-L/2, L/2, M)';
    y0 = linspace(-L/2, L/2, M)';
    dx = x0(2) - x0(1);
    dy = y0(2) - y0(1);
    dkx = 2 * pi / M / dx;
    dky = 2 * pi / M / dy;
    kx0 = (-M/2:M/2-1)' * dkx;
    ky0 = (-M/2:M/2-1)' * dky;
    [meshx, meshy] = meshgrid(x0, y0);
    [meshkx, meshky] = meshgrid(kx0, ky0);
    % construct TU on k grid
    T = (meshkx.^2 + meshky.^2) / 2 / mass;
    TU = exp(-1i * dt * T);
    TU = fftshift(TU);
    % construct VU
    VU = zeros(M,M,N,N);
    Hs = zeros(M,M,N,N);
    for j=1:M
        for k=1:M
            A = 0.01;
            B = 1.6;
            C = 0.005;
            D = 1.0;
            W = 0.0;
            phi =  W * y0(k);
            eip = exp(1i * phi);

            H = zeros(2,2);
            x = x0(j);
            if x >= 0.0
                H(1,1) = A * (1.0 - exp(-B * x));
            else
                H(1,1) = -A * (1.0 - exp(B * x));
            end
            H(2,2) = -H(1,1);
            H(1,2) = C * exp(-D * x * x) * eip;
            H(2,1) = conj(H(1,2));
            VU(k,j,:,:) = expm(-1i * dt / 2 * H);
            Hs(k,j,:,:) = H;
        end
    end
    % construct phase corrected eva & evt, optional
    %{
    evas = zeros(M,M,N,N);
    evts = zeros(M,M,N,N);
    for j=1:M
        for k=1:M
            [evt, eva] = eig(reshape(Hs(k,j,:,:), 2, 2));
            if j == 1 && k == 1
                evts(k,j,:,:) = evt;
                evas(k,j,:,:) = eva;
            elseif k == 1
                lastevt = reshape(evts(k,j-1,:,:), 2, 2);
                phase = lastevt' * evt;
                evt(:,1) = evt(:,1) / phase(1,1);
                evt(:,2) = evt(:,2) / phase(2,2);
                evts(k,j,:,:) = evt;
                evas(k,j,:,:) = eva;
            else
                lastevt = reshape(evts(k-1,j,:,:), 2, 2);
                phase = lastevt' * evt;
                evt(:,1) = evt(:,1) / phase(1,1);
                evt(:,2) = evt(:,2) / phase(2,2);
                evts(k,j,:,:) = evt;
                evas(k,j,:,:) = eva;
            end
        end
    end
    %}
    % Initial wavefunction -- Gaussian wavepacket
    psi0 = zeros(M,M,N);
    psi0(:,:,1) = c1 * exp(1i*(kxI*meshx + kyI*meshy)) .* exp(-(meshx-xI).^2/sigmax^2 - (meshy-yI).^2/sigmay^2);
    psi0(:,:,2) = c2 * exp(1i*(kxI*meshx + kyI*meshy)) .* exp(-(meshx-xI).^2/sigmax^2 - (meshy-yI).^2/sigmay^2);
    psi0 = psi0 / sqrt(sum(sum(sum(abs(psi0).^2))));
    % psim & psi_k_m -- for plot
    psi0_k(:,:,1) = fftshift(fft2(psi0(:,:,1)));
    psi0_k(:,:,2) = fftshift(fft2(psi0(:,:,2)));
    psi0_k = psi0_k / sqrt(sum(sum(sum(abs(psi0_k).^2))));
    psi_max = max(max(max(abs(psi0).^2)));
    psi_k_max = max(max(max(abs(psi0_k).^2)));
    % propagate WF
    ana_step = 0;
    psi = psi0;
    for t=0:Nstep-1
        % exp(-iVdt/2) * |Psi> in diab repr
        psi_k(:,:,1) = VU(:,:,1,1).*psi(:,:,1) + VU(:,:,1,2).*psi(:,:,2);
        psi_k(:,:,2) = VU(:,:,2,1).*psi(:,:,1) + VU(:,:,2,2).*psi(:,:,2);
        % exp(-iTdt) * psi
        psi_k(:,:,1) = TU .* fft2(psi_k(:,:,1));
        psi_k(:,:,2) = TU .* fft2(psi_k(:,:,2));
        % exp(-iVdt/2) * psi
        psi_k(:,:,1) = ifft2(psi_k(:,:,1));
        psi_k(:,:,2) = ifft2(psi_k(:,:,2));
        psi(:,:,1) = VU(:,:,1,1).*psi_k(:,:,1) + VU(:,:,1,2).*psi_k(:,:,2);
        psi(:,:,2) = VU(:,:,2,1).*psi_k(:,:,1) + VU(:,:,2,2).*psi_k(:,:,2);
        % analysis & report
        if mod(t,tgraph) == 0
            % analysis in diab
            psi_k(:,:,1) = fftshift(fft2(psi(:,:,1)));
            psi_k(:,:,2) = fftshift(fft2(psi(:,:,2)));
            psi_k = psi_k / sqrt(sum(sum(sum(abs(psi_k).^2))));
            
            norm_k1 = sum(sum(abs(psi_k(:,:,1)).^2)) + 1e-16;
            norm_k2 = sum(sum(abs(psi_k(:,:,2)).^2)) + 1e-16;

            p1x = sum(sum(abs(psi_k(:,:,1)).^2 .* meshkx)) / norm_k1;
            p1y = sum(sum(abs(psi_k(:,:,1)).^2 .* meshky)) / norm_k1;
            p2x = sum(sum(abs(psi_k(:,:,2)).^2 .* meshkx)) / norm_k2;
            p2y = sum(sum(abs(psi_k(:,:,2)).^2 .* meshky)) / norm_k2;

            KE = sum(sum((abs(psi_k(:,:,1)).^2 + abs(psi_k(:,:,2)).^2) .* (meshkx.^2 + meshky.^2) / 2 / mass));
            PE = sum(sum( conj(psi(:,:,1)) .* Hs(:,:,1,1) .* psi(:,:,1) + conj(psi(:,:,1)) .* Hs(:,:,1,2) .* psi(:,:,2) + conj(psi(:,:,2)) .* Hs(:,:,2,1) .* psi(:,:,1) + conj(psi(:,:,2)) .* Hs(:,:,2,2) .* psi(:,:,2) ));

            % analysis in adiab
            %{
            psiad(:,:,1) = conj(evts(:,:,1,1)).*psi(:,:,1)+conj(evts(:,:,2,1)).*psi(:,:,2);
            psiad(:,:,2) = conj(evts(:,:,1,2)).*psi(:,:,1)+conj(evts(:,:,2,2)).*psi(:,:,2);
            psiad = psiad / sqrt(sum(sum(sum(abs(psiad).^2))));
            psiad_k(:,:,1) = fftshift(fft2(psiad(:,:,1)));
            psiad_k(:,:,2) = fftshift(fft2(psiad(:,:,2)));
            psiad_k = psiad_k / sqrt(sum(sum(sum(abs(psiad_k).^2))));

            normad_k1 = sum(sum(abs(psiad_k(:,:,1)).^2)) + 1e-16;
            normad_k2 = sum(sum(abs(psiad_k(:,:,2)).^2)) + 1e-16;

            pad1x = sum(sum(abs(psiad_k(:,:,1)).^2 .* meshkx)) / normad_k1;
            pad1y = sum(sum(abs(psiad_k(:,:,1)).^2 .* meshky)) / normad_k1;
            pad2x = sum(sum(abs(psiad_k(:,:,2)).^2 .* meshkx)) / normad_k2;
            pad2y = sum(sum(abs(psiad_k(:,:,2)).^2 .* meshky)) / normad_k2;
            %}

            % output
            if t == 0
                fprintf(fID, '# EXACT\n');
                fprintf(fID, '# xI = %8.4f yI = %8.4f kxI = %8.4f kyI = %8.4f sigmax = %8.4f sigmay = %8.4f A = %8.4f init_s = %8.4f c1 = %8.4f c2 = %8.4f \n', ...
                                    xI, yI, kxI, kyI, sigmax, sigmay, A, init_s, c1, c2);
                fprintf(fID, '# L = %8.4f M = %8d dt = %8.4f Nstep = %8d tgraph = %8d\n', ...
                                    L, M, dt, Nstep, tgraph);
                fprintf(fID, '#%9s%16s%16s%16s%16s%16s%16s%16s\n', ...
                                't', 'n0', 'n1', 'px0', 'px1', 'py0', 'py1', 'Etot');
            end

            %{
            fprintf(fID, '#%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f\n', ...
                        t*dt, ...
                        sum(sum(sum(abs(psiad(:,:,1).^2)))), ...
                        sum(sum(sum(abs(psiad(:,:,2).^2)))), ...
                        pad1x, pad2x, pad1y, pad2y, ...
                        KE + PE ...
                        );
            %}

            fprintf(fID, '#%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f\n', ...
                        t*dt, ...
                        sum(sum(sum(abs(psi(:,:,1).^2)))), ...
                        sum(sum(sum(abs(psi(:,:,2).^2)))), ...
                        p1x, p2x, p1y, p2y, ...
                        KE+PE ...
                        );
            % plot
            if enable_plot == true
                subplot(2,2,1);
                contour(y0,x0,abs(psi(:,:,1)).^2,[psi_max/100:psi_max/100:psi_max]);
                title('Real space -- Pop Diab 1');

                subplot(2,2,2);
                contour(y0,x0,abs(psi(:,:,2)).^2,[psi_max/100:psi_max/100:psi_max]);
                title('Real space -- Pop Diab 2');

                subplot(2,2,3);
                contour(ky0,kx0,abs(psi_k(:,:,1)).^2,[psi_k_max/100:psi_k_max/100:psi_k_max]);
                title('Mom space -- Pop Diab 1');

                subplot(2,2,4);
                contour(ky0,kx0,abs(psi_k(:,:,2)).^2,[psi_k_max/100:psi_k_max/100:psi_k_max]); 
                title('Mom space -- Pop Diab 2');
                drawnow;
            end

        end
    end
    fprintf(fID, '%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f%16.10f\n', ...
            kxI, ...
            sum(sum(sum(abs(psi(:,:,1).^2)))), ...
            sum(sum(sum(abs(psi(:,:,2).^2)))), ...
            p1x, p2x, p1y, p2y ...
            );
    fprintf('# '); toc;
end
