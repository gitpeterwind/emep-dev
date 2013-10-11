! DO3SE assertion helpers, depends on DO3SE_Util_ml
#define ASSERT(check) call assert((check), "ASSERTION FAILED: check")
#define ASSERT_DEFINED(id) call assert(is_def(id), "id not defined")
#define UNKNOWN_STRING(id) call assert(.false., "unknown id: "//trim((id)))
#define ERROR(msg) call assert(.false., (msg))
