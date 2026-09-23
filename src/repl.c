#define _POSIX_C_SOURCE 200809L

#include "smallcas.h"

#include <readline/history.h>
#include <readline/readline.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int sc_is_exit_command(const char *line)
{
    return strcmp(line, "quit") == 0 || strcmp(line, "exit") == 0;
}

static int sc_is_blank(const char *line)
{
    while (*line == ' ' || *line == '\t')
        line++;
    return *line == '\0';
}

static char *sc_read_line(int interactive)
{
    char *line = NULL;
    size_t cap = 0;
    ssize_t n;

    if (interactive)
        return readline("smallcas> ");
    n = getline(&line, &cap, stdin);
    if (n < 0) {
        free(line);
        return NULL;
    }
    if (n != 0 && line[n - 1] == '\n')
        line[n - 1] = '\0';
    return line;
}

static void sc_history_filename(char *path, size_t size)
{
    const char *name = getenv("SMALLCAS_HISTORY");
    const char *home = getenv("HOME");

    if (name != NULL && name[0] != '\0')
        snprintf(path, size, "%s", name);
    else if (home != NULL && home[0] != '\0')
        snprintf(path, size, "%s/.smallcas_history", home);
    else
        path[0] = '\0';
}

static void sc_history_start(const char *path)
{
    rl_readline_name = "smallcas";
    using_history();
    stifle_history(10000);
    if (path[0] != '\0')
        read_history(path);
}

static void sc_history_add(const char *path, const char *line)
{
    add_history(line);
    if (path[0] != '\0')
        write_history(path);
}

int main(void)
{
    sc_context ctx;
    char history_path[4096];
    int interactive = isatty(STDIN_FILENO);
    char *line;

    sc_context_init(&ctx);
    sc_history_filename(history_path, sizeof(history_path));
    if (interactive)
        sc_history_start(history_path);
    puts("smallcas iter30 -- ZZ and univariate ZZ polynomials -- type quit to exit");
    while ((line = sc_read_line(interactive)) != NULL) {
        if (sc_is_blank(line)) {
            free(line);
            continue;
        }
        if (interactive)
            sc_history_add(history_path, line);
        if (sc_is_exit_command(line)) {
            free(line);
            break;
        }
        if (!sc_repl_line(&ctx, line))
            fprintf(stderr, "error: %s\n", ctx.error);
        else if (ctx.show_result)
            sc_value_print(ctx.result);
        free(line);
    }
    rl_clear_history();
    sc_context_clear(&ctx);
    return 0;
}
