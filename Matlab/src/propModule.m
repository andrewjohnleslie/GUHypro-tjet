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
classdef propModule < handle
    properties (SetAccess=protected)
        this = 0
        N1_
        N2_
        owner=true
    end
    methods 
        function obj = propModule(in,varargin)
           if nargin>0
                 if isa(in,'uint64')
                    obj.this = in;
                    if nargin>1
                        if ~islogical(varargin{1})
                            error('Error: the second input must be a boolean');
                        end
                        obj.owner = varargin{1};
                    end
                else
                    error('Input must be either a mixture or a uint32');
                end
            end
        end
        function delete(obj)
            if obj.this~=0
                HyPro.propModuleMex('destructor',obj.this);
                obj.this = 0;
                
                obj.N1_.delete;
                obj.N2_.delete;
            end
        end
        function out = Name(obj,varargin)
            if nargin==1
                out = cell(size(obj));
                for i=1:length(obj)
                    out{i} = HyPro.propModuleMex('getName',obj(i).this);
                end
                if length(obj)==1
                    out = out{1};
                end
            else
                HyPro.propModuleMex('setName',obj.this,varargin{1});
            end
        end
        function verbosity(obj,verb)
            HyPro.propModuleMex('setVerbosity',obj.this,verb);
        end
        function out = calculate(obj)
            out = HyPro.propModuleMex('calculate',obj.this);
        end
        function unreduce(obj)
            HyPro.propModuleMex('unreduce',obj.this);
        end
    end
    methods(Access=protected)
        function N1(obj)
            obj.N1_ = HyPro.Node(HyPro.propModuleMex('getN1',obj.this));
        end
        function N2(obj)
            obj.N2_ = HyPro.Node(HyPro.propModuleMex('getN2',obj.this));
        end
    end
end
