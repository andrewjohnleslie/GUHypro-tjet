/*!
 *  \author     Alessandro Mogavero
 *  \n
 *  \email      alessandro.mogavero@strath.ac.uk
 *  \version    1.0
 *  \copyright  Copyright 2016 Alessandro Mogavero
 */
 /* This file is part of HyPro.
 *
 * HyPro is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HyPro is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HyPro.  If not, see <http://www.gnu.org/licenses/>. */
classdef systemModule < HyPro.propModule
    methods 
        function obj = systemModule(fileName,specieGas)
%             prevDir = cd;
%             cd /home/trb12187/Matlab/SpaceplaneModel/HyPro/Test/GP
%             obj.this = GP(varargin{:});
%             cd(prevDir);

%             [~,~,~,obj.this] = Hyperion(1e5,300,2);
            if ~ischar(fileName)
                error('HYPRO:invalidIO','Error: fileName must be a string');
            end
            if ~isa(specieGas,'HyPro.mixture')
                error('HYPRO:invalidIO','Error: specieGas must be a HyPro.mixture.');
            end
            obj.this = HyPro.systemModuleMex('constructor',fileName,specieGas.this);

            obj.N1;
            obj.N2;
        end
        function out = nodes(obj)
            h = HyPro.systemModuleMex('nodes',obj.this);
            for i=1:length(h)
                out(i) = HyPro.Node(h(i),false);
            end
        end
        function out = thrust(obj)
            out = HyPro.systemModuleMex('thrust',obj.this);
        end
        function out = propellantMfr(obj)
            out = HyPro.systemModuleMex('propellantMfr',obj.this);
        end
        function out = modules(obj)
            h = HyPro.systemModuleMex('modules',obj.this);
            for i=1:length(h)
                out(i) = HyPro.propModule(h(i),false);
            end
        end
    end
end
