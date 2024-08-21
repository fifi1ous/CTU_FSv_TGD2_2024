function [Nav_rin] = export_cell(u)

%                      1            2      3       4       5       6        7               8       9        10      11       12          13       14       15      16             
Nav_rin=[u{5},u{2},u{8}(1)/1000,u{8}(2),u{8}(3),u{8}(4),u{8}(5),u{8}(6),(u{8}(7)^2)/1000,u{8}(8),u{8}(9),u{8}(10),u{8}(11),u{8}(12),u{8}(13)/1000,u{8}(14),u{8}(15),u{8}(16)];
end