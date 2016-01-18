help_lang_dir = get_absolute_file_path("build_help.sce");
TOOLBOX_TITLE = "fgoalattain"
tbx_build_help(TOOLBOX_TITLE, help_lang_dir);
ok = add_help_chapter("Demo",get_absolute_file_path("build_help.sce"));
clear help_lang_dir;
