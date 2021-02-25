classdef LutModel
    %LutModel lookuptable for complex numbers
    properties
        xin
        yout
        interptype
        casc
        noise
    end
    
    methods
        function obj = LutModel(amplin, iqout, interptype)
            %LUTMODEL Construct an instance of this class.
            obj.xin = amplin; %Amplitude column of the LUT
            obj.yout = iqout;
            obj.casc = [];
            obj.noise = 0;
            if nargin == 3
                obj.interptype = interptype;
            else
                obj.interptype = 'spline';
            end
        end
        function pa_output = transmit(obj, in)
            %transmit Send through the LUT
            %         The LUT is amplitude in, I, Q.
            %
            yint = interp1(obj.xin, obj.yout, abs(in), obj.interptype, obj.yout(end));
            %exp(1j*phase(in)); %this takes too long
            normalized = in ./ abs(in);
            pa_output = normalized.*yint + rand(length(yint),1)*obj.noise;
            
            % do the cascading. Other objects can be cascaded, keep in
            % mind, and there can be multiple objects in a single cascade.
            % this is a recipe for infinite loops of course.
            for i = 1:length(obj.casc)
                pa_output = obj.casc(i).transmit(pa_output);
            end
        end
        
        function inv = invert(obj, npoints)
            %invert Invert the output of the LUT
            %       Can be used to make a predistorter or postdistorter.
            % Scale the input to correspond to the output amplitude.
            % This means that we expect the output max amplitude to be the
            % same as the input max amplitude.
            % We should try to figure out a better way to scale data, then.
            % That said if you keep this method, you're guaranteed never
            % to clip: louder data will just get clipped to the peak LUT
            % ampl.
            outampls = abs(obj.yout);

            if nargin == 2
                tablesize = npoints;
            else
                [tablesize,~] = size(outampls);
            end

            new_x_basis = linspace(min(outampls), max(outampls), tablesize);
            new_yout = interp1(outampls, obj.xin, new_x_basis, obj.interptype);

            inv = LutModel(new_x_basis, new_yout, obj.interptype);
            
            for i = 1:length(obj.casc)
                casc_inv = obj.casc(i).invert();
                inv = casc_inv.cascade(inv);
            end
        end
        
        function obj = cascade(obj, tbc)
            %cascade Cascade the LUT with another one.
            % This is used to simulate transmission through two or more
            % consecutive LUTs.
            
            %we can just pass the output of one into the other
            %the arg is the second one.
            obj.casc = horzcat(obj.casc, tbc);
        end
    end
            
end
