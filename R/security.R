#' Get Machine Code
#'
#' @return Formatted machine identifier (xxxx-xxxx-xxxx-xxxx)
#' @export
get_machine_code <- function() {
    # Combine system information to generate a unique identifier
    sys_info <- Sys.info()
    # Attempt to get a more stable hardware identifier
    hw_id <- tryCatch(
        {
            if (Sys.info()["sysname"] == "Darwin") {
                # macOS
                res <- system("ioreg -rd1 -c IOPlatformExpertDevice | grep -E 'IOPlatformUUID'", intern = TRUE)
                gsub(".*\"IOPlatformUUID\" = \"(.*)\"", "\\1", res)
            } else if (Sys.info()["sysname"] == "Linux") {
                # Linux
                if (file.exists("/etc/machine-id")) {
                    readLines("/etc/machine-id", n = 1)
                } else {
                    paste(sys_info, collapse = "")
                }
            } else {
                # Windows
                res <- system("wmic csproduct get uuid", intern = TRUE)
                paste(res, collapse = "")
            }
        },
        error = function(e) paste(sys_info, collapse = "")
    )

    raw_str <- paste0(hw_id, sys_info["user"], sys_info["nodename"])
    # MD5 Digest
    hash <- digest::digest(raw_str, algo = "md5")

    # Format into xxxx-xxxx-xxxx-xxxx (32 chars md5 -> pick 16 or use all?)
    # Let's take groups of 4 from the hash to make it look like a key
    formatted <- paste(
        substr(hash, 1, 4),
        substr(hash, 5, Group2 <- 8),
        substr(hash, 9, 12),
        substr(hash, 13, 16),
        sep = "-"
    )
    return(toupper(formatted))
}

#' Verify Activation Code
#'
#' @param activation_code Activation code entered by user
#' @return Logical value
#' @keywords internal
verify_activation <- function(activation_code) {
    if (missing(activation_code) || is.null(activation_code) || activation_code == "") {
        return(FALSE)
    }

    machine_code <- get_machine_code()
    # MD5(machine_code_formatted + salt)
    salt <- "SMRR_SECURE_SALT_2026_DTN_V2"
    expected_code <- digest::digest(paste0(machine_code, salt), algo = "md5")
    # Format activation code too for consistency if needed, but let's keep it simple

    return(activation_code == expected_code)
}

#' Get License File Path
#' @keywords internal
get_license_path <- function() {
    file.path(Sys.getenv("HOME"), ".smrr_license")
}

#' Activate SMRr Package
#'
#' @description
#' Enter the activation code provided by the author to permanently activate the current machine.
#'
#' @param activation_code Activation code string
#' @return TRUE on success, otherwise stops with error
#' @export
activate_smrr <- function(activation_code) {
    if (verify_activation(activation_code)) {
        writeLines(activation_code, get_license_path())
        message("激活成功！您的 SMRr 包已永久激活。")
        return(invisible(TRUE))
    } else {
        stop(sprintf("\n[激活失败 / Activation Failed]\n无效的激活码！\n您的机器码是：%s\n\n请联系作者获取激活码：\n邮箱：202201230726@163.com\n微信：CangMing-03", get_machine_code()))
    }
}

#' Check Activation Status
#' @param activation_code Optional activation code for one-time use
#' @keywords internal
check_auth <- function(activation_code = NULL) {
    # 1. Check provided activation_code
    if (!is.null(activation_code) && verify_activation(activation_code)) {
        return(TRUE)
    }

    # 2. Check local license file
    lic_path <- get_license_path()
    if (file.exists(lic_path)) {
        saved_code <- readLines(lic_path, n = 1, warn = FALSE)
        if (verify_activation(saved_code)) {
            return(TRUE)
        }
    }

    # 3. Verification failed
    stop(sprintf("\n[软件未激活 / SMRr Not Activated]\n请先运行：activate_smrr('您的激活码') 进行永久激活。\n您的机器码是：%s\n\n联系作者获取激活码：\n邮箱：202201230726@163.com\n微信：CangMing-03", get_machine_code()))
}
