function [xhat, meas] = gyro_mag_no_acc(calAcc, calGyr, calMag)
% FILTERTEMPLATE  Filter template
%
% This is a template function for how to collect and filter data
% sent from a smartphone live.  Calibration data for the
% accelerometer, gyroscope and magnetometer assumed available as
% structs with fields m (mean) and R (variance).
%
% The function returns xhat as an array of structs comprising t
% (timestamp), x (state), and P (state covariance) for each
% timestamp, and meas an array of structs comprising t (timestamp),
% acc (accelerometer measurements), gyr (gyroscope measurements),
% mag (magnetometer measurements), and orint (orientation quaternions
% from the phone).  Measurements not availabe are marked with NaNs.
%
% As you implement your own orientation estimate, it will be
% visualized in a simple illustration.  If the orientation estimate
% is checked in the Sensor Fusion app, it will be displayed in a
% separate view.
%
% Note that it is not necessary to provide inputs (calAcc, calGyr, calMag).

  %% Setup necessary infrastructure
  import('com.liu.sensordata.*');  % Used to receive data.

  %% Filter settings
  t0 = [];  % Initial time (initialize on first data received)
  nx = 4;   % Assuming that you use q as state variable.
  % Add your filter settings here.
  R_acc = [4.05693372065651e-05	-3.99459058972891e-07	-2.78969000314148e-06
-3.99459058972891e-07	4.09687938838409e-05	1.85491365223252e-06
-2.78969000314148e-06	1.85491365223252e-06	6.26278673591397e-05];
  R_gyr = [5.46894781840160e-07	-2.28728996080674e-09	-1.08637456955466e-07
-2.28728996080674e-09	4.12495719075584e-07	4.21328806218598e-08
-1.08637456955466e-07	4.21328806218598e-08	3.89249867409969e-07];
  R_mag = [0.0953350666994726	-0.00334500918359894	0.00947364856139815
-0.00334500918359894	0.102573124990338	-0.00538472040872685
0.00947364856139815	-0.00538472040872685	0.109661863216162];
  g0 = [0.0665 0.0097 9.9057]';
  acc_range = 0.15;%0.8acc-1.2acc is available
  alpha = 0.01;%alpha for mag
  mag_range = 0.10;
% 2.97686651183955
% 18.4024368960637
% -40.5299425720934
  m0 = [0 18.6417 -40.5299]';
  L = norm(m0);
  % Current filter state.
  x = [1; 0; 0 ;0];
  P = eye(nx, nx);

  % Saved filter states.
  xhat = struct('t', zeros(1, 0),...
                'x', zeros(nx, 0),...
                'P', zeros(nx, nx, 0));

  meas = struct('t', zeros(1, 0),...
                'acc', zeros(3, 0),...
                'gyr', zeros(3, 0),...
                'mag', zeros(3, 0),...
                'orient', zeros(4, 0));
  try
    %% Create data link
    server = StreamSensorDataReader(3400);
    % Makes sure to resources are returned.
    sentinel = onCleanup(@() server.stop());

    server.start();  % Start data reception.

    % Used for visualization.
    figure(1);
    subplot(1, 2, 1);
    ownView = OrientationView('Own filter', gca);  % Used for visualization.
    googleView = [];
    counter = 0;  % Used to throttle the displayed frame rate.

    %% Filter loop
    while server.status()  % Repeat while data is available
      % Get the next measurement set, assume all measurements
      % within the next 5 ms are concurrent (suitable for sampling
      % in 100Hz).
      data = server.getNext(5);

      if isnan(data(1))  % No new data received
        continue;        % Skips the rest of the look
      end
      t = data(1)/1000;  % Extract current time

      if isempty(t0)  % Initialize t0
        t0 = t;
      end
      acc = data(1, 2:4)';
      if ~any(isnan(acc))  % Acc measurements are available.
%           if abs(norm(acc)-norm(g0)) < norm(g0)*acc_range
%               [x, P] = mu_g(x,P,acc,R_acc,g0);
%               [x, P] = mu_normalizeQ(x, P);
%               ownView.setAccDist(0);
%           else
%               ownView.setAccDist(1);
%           end

      end
     
      gyr = data(1, 5:7)';
      if ~any(isnan(gyr))  % Gyro measurements are available.
        [x, P] = tu_qw(x,P,gyr,0.001, R_gyr);
        [x, P] = mu_normalizeQ(x,P);
      else
        [x, P] = tu_qw_no(x,P,gyr,0.001, R_gyr);
        [x, P] = mu_normalizeQ(x,P);
      end


      mag = data(1, 8:10)';
      if ~any(isnan(mag))  % Mag measurements are available.
        L = (1-alpha)*L + alpha*norm(mag);
        if abs(L-norm(mag)) < L*mag_range
           [x,P] = mu_m(x,P,mag,m0,R_mag);
           [x,P] = mu_normalizeQ(x,P);
           ownView.setMagDist(0);
        else
           ownView.setMagDist(1);
        end
      end

      orientation = data(1, 18:21)';  % Google's orientation estimate.

      % Visualize result
      if rem(counter, 10) == 0
        setOrientation(ownView, x(1:4));
        title(ownView, 'OWN', 'FontSize', 16);
        if ~any(isnan(orientation))
          if isempty(googleView)
            subplot(1, 2, 2);
            % Used for visualization.
            googleView = OrientationView('Google filter', gca);
          end
          setOrientation(googleView, orientation);
          title(googleView, 'GOOGLE', 'FontSize', 16);
        end
      end
      counter = counter + 1;

      % Save estimates
      xhat.x(:, end+1) = x;
      xhat.P(:, :, end+1) = P;
      xhat.t(end+1) = t - t0;

      meas.t(end+1) = t - t0;
      meas.acc(:, end+1) = acc;
      meas.gyr(:, end+1) = gyr;
      meas.mag(:, end+1) = mag;
      meas.orient(:, end+1) = orientation;
    end
  catch e
    fprintf(['Unsuccessful connecting to client!\n' ...
      'Make sure to start streaming from the phone *after*'...
             'running this function!']);
  end
end
