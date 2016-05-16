
package edu.usc.scrc.PriorityPruner;

import java.security.Permission;

/***
 * From http://neverfear.org/blog/view/157/Testing_code_that_calls_System_exit_in_Java
 *
 */
public class ExitDeniedSecurityManager extends SecurityManager {
	 
    public static final class ExitSecurityException extends SecurityException {
        private final int status;
 
        public ExitSecurityException(final int status) {
            this.status = status;
        }
 
        public int getStatus() {
            return this.status;
        }
    }
 
    @Override
    public void checkExit(final int status) {
        throw new ExitSecurityException(status);
    }
 
    @Override
    public void checkPermission(final Permission perm) {}
 
    /* TODO: Make sure you override everything with an empty method as above */
}