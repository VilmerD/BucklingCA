%% Set up cluster
c = parcluster;
c.AdditionalProperties.ProjectName = 'lu2022-2-24';
c.AdditionalProperties.WallTime = '00:10:00';
c.saveProfile;