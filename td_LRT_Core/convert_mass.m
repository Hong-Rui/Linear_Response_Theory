function [mass] = convert_mass(mass_text_filename)
	fileID = fopen(mass_text_filename,'r');
	mass = textscan(fileID,'%.12f');
	fclose(fileID);

	mass = cell2mat(mass);

end
