clc; clearvars; close all 
% Carregar os arquivos
refFile = load('montessoriphistorage.mat'); % Arquivo de referência
testFile = load('testphistorage.mat');      % Arquivo de teste 1
vectorFile = load('vectorphistorage.mat');  % Arquivo de teste 2

% Extraia os arrays
phi_ref = refFile.phistorage;
phi_test = testFile.phistorage;
phi_vector = vectorFile.phistorage;

% Verificar se os arrays são idênticos
are_test_identical = isequal(phi_ref, phi_test);
are_vector_identical = isequal(phi_ref, phi_vector);

% Exibir os resultados
if are_test_identical
    disp('testphistorage.mat é idêntico a montessoriphistorage.mat.');
else
    disp('testphistorage.mat NÃO é idêntico a montessoriphistorage.mat.');
end

if are_vector_identical
    disp('vectorphistorage.mat é idêntico a montessoriphistorage.mat.');
else
    disp('vectorphistorage.mat NÃO é idêntico a montessoriphistorage.mat.');
end
